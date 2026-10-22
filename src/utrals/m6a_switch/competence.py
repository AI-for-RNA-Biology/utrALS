#!/usr/bin/env python
"""
For each proximal OR distal isoform in table_of_pairs.csv (--isoform):
  1. Locate its anchor in the 3'UTR GTF annotation (--anchor): the PAS (3' end of the 3'UTR
     isoform, default) or the 3'UTR start (5' end of the record = first nt after the stop codon).
  2. Extract the [WL:WR) nt window around that anchor (5'->3' transcript orientation; --wl/--wr,
     e.g. 0 200 = the 200 nt downstream of the proximal PAS, -350 0 = the 350 nt upstream of
     the distal PAS, --anchor utr_start -100 100 = 100 nt of CDS (incl. the stop codon) + the
     first 100 nt of the 3'UTR). Position 0 = the first nt after the PAS / the first 3'UTR nt.
     With --cut-at-exon-start (utr_start only), a window whose upstream part would leave the exon
     containing the stop codon (exon structure from GENCODE v44, matched by stop-codon position)
     starts at that exon's first nt instead, so no intron is folded; the per-transcript window
     start is saved as per_transcript.csv window_wl and reused by every downstream step.
  3. Fold the window with the ViennaRNA partition function to get base-pair probabilities.
  4. For every measured m6A site (RMBase) landing in that window, find its most probable
     pairing partner j and test whether it is "switch-competent":
       p(i, j) >= 0.5  AND  j itself is part of a run of >= 3 consecutive U's (no flank
       window around j -- j must directly be one of the U's in that run).
  5. Aggregate switch-competent site counts/fractions per transcript, with the
     nuclear-lengthened label of the pair (pairs_distals_D22_Nuclear_Wildtype_vs_D35_Nuclear_Wildtype.csv).

Run with no arguments to reproduce the original proximal [0:200] run (legacy defaults, legacy
output directory m6A_switch_competence_0_200). Every run writes params.json next to its
outputs (arguments, input paths, script checksums, package versions) for provenance.

Outputs under --out-dir:
  - per_site.csv       : one row per (transcript, m6A site) found in the window. Key column is
                          "<isoform>_id". u_tract_genomic_start/end are 0-based and INCLUSIVE.
  - per_transcript.csv : one row per transcript (all isoforms, including those with no m6A
                          site), with switch-competent counts/fractions, the window's A count
                          (n_A, for A-normalized statistics) and the nuclear-lengthened label
  - params.json
"""

import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
import RNA

from ..paths import ANNOT_DIR, REPO_ROOT, PROJECT_DIR as _PROJECT_DIR  # data locations, see paths.env

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TABLE_OF_PAIRS = f"{_PROJECT_DIR}/MN_diff/results/3UTR_analysis/table_of_pairs.csv"
SHIFT_FILE = f"{_PROJECT_DIR}/MN_diff/results/3UTR_analysis/shifts/pairs_distals_D22_Nuclear_Wildtype_vs_D35_Nuclear_Wildtype.csv"
GTF_PATH = f"{ANNOT_DIR}/3UTR_MNs/hg38_Denovo_Ensembl_RefSeq_Filt_Reformat.gtf"
M6A_BED = f"{ANNOT_DIR}/RMBase/human.hg38.m6A.result.col29.bed"
GENOME_FASTA = f"{ANNOT_DIR}/GENCODE/GRCh38.p14.genome.fa"
GENCODE_GTF = f"{ANNOT_DIR}/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"

M6A_RESULTS_ROOT = f"{_PROJECT_DIR}/MN_diff/results/RNA_modifs/m6A"

# Legacy defaults (proximal isoform, [0:200) downstream of its PAS). Kept as module-level names
# because other scripts import them -- but new code should pass wl/wr explicitly instead of
# relying on these, so a change here can never silently change another analysis' window.
WL, WR = 0, 200
WIN_LEN = WR - WL  # 200

OUT_DIR = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_{WL}_{WR}"

P_THRESHOLD = 0.5
U_TRACT_MIN = 3


# ---------------------------------------------------------------------------
# GTF -> one row per transcript (chrom, start, end, strand), PAS position
# ---------------------------------------------------------------------------
def load_pas_lookup(gtf_path, needed_ids=None):
    """Parse the GTF into one row per transcript (chrom, start, end, strand, PAS).
    If `needed_ids` is given, pre-filter the file with `grep -F` (fixed-string,
    one pattern per transcript id) so we never parse/allocate rows for the
    ~2.5M transcripts in the full annotation that we don't need -- table_of_pairs
    only references a few thousand of them."""
    rows = []
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    gene_re = re.compile(r'gene_id "([^"]+)"')

    if needed_ids:
        import subprocess
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            tf.write("\n".join(needed_ids))
            pattern_file = tf.name
        try:
            proc = subprocess.run(
                ["grep", "-F", "-f", pattern_file, gtf_path],
                capture_output=True, text=True, check=False,
            )
        finally:
            os.remove(pattern_file)
        lines = proc.stdout.splitlines()
    else:
        with open(gtf_path) as f:
            lines = f.readlines()

    needed_ids = set(needed_ids) if needed_ids else None
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        chrom, _, _, start, end, _, strand, _, attrs = fields
        tx_match = tx_re.search(attrs)
        if not tx_match:
            continue
        transcript_id = tx_match.group(1)
        if needed_ids is not None and transcript_id not in needed_ids:
            continue
        rows.append((
            gene_re.search(attrs).group(1) if gene_re.search(attrs) else None,
            transcript_id,
            chrom if chrom.startswith("chr") else "chr" + chrom,
            int(start), int(end), strand,
        ))
    df = pd.DataFrame(rows, columns=["gene_id", "transcript_id", "chrom", "start", "end", "strand"])
    df = df.drop_duplicates(subset="transcript_id", keep="first")
    df["PAS"] = np.where(df["strand"] == "+", df["end"], df["start"])
    # 3'UTR start in the same convention as PAS: window [wl, wr) -> + strand [a+wl, a+wr),
    # - strand [a-wr, a-wl); position 0 = first 3'UTR nt (0-based start-1 on +, end-1 on -)
    df["UTR_START"] = np.where(df["strand"] == "+", df["start"] - 1, df["end"])
    return df.set_index("transcript_id")[["gene_id", "chrom", "strand", "PAS", "UTR_START"]]


ANCHORS = ("pas", "utr_start")


def anchor_coordinate(info, anchor):
    """Window anchor of one transcript (row of load_pas_lookup) for --anchor pas / utr_start."""
    return int(info["PAS"] if anchor == "pas" else info["UTR_START"])


def stop_exon_nt_before_utr(lookup, gtf_path=GENCODE_GTF):
    """For each transcript of `lookup` (load_pas_lookup output), the number of nt of the exon
    containing its stop codon that lie upstream of the first 3'UTR nt (stop codon included),
    from GENCODE: transcripts are matched on the stop codon ending right before the 3'UTR start;
    when several GENCODE transcripts share it, the longest such exon part is used. NaN if no match."""
    chroms = set(lookup["chrom"])
    stops, exons = [], []
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t", 8)
            if x[2] not in ("stop_codon", "exon") or x[0] not in chroms:
                continue
            row = (tid_re.search(x[8]).group(1), x[0], int(x[3]), int(x[4]), x[6])
            (stops if x[2] == "stop_codon" else exons).append(row)
    cols = ["tid", "chrom", "s", "e", "strand"]
    st, ex = pd.DataFrame(stops, columns=cols), pd.DataFrame(exons, columns=cols)
    st["utr0"] = np.where(st["strand"] == "+", st["e"], st["s"] - 2)          # 0-based first 3'UTR nt
    ex = ex.merge(st[["tid", "s", "e", "utr0"]].rename(columns={"s": "ss", "e": "se"}), on="tid")
    ex = ex[(ex["s"] <= ex["ss"]) & (ex["e"] >= ex["se"])]                   # exon holding the stop codon
    ex["nt_before_utr"] = np.where(ex["strand"] == "+", ex["se"] - ex["s"] + 1, ex["e"] - ex["ss"] + 1)
    best = ex.groupby(["chrom", "strand", "utr0"])["nt_before_utr"].max()
    utr0 = np.where(lookup["strand"] == "+", lookup["UTR_START"], lookup["UTR_START"] - 1)
    keys = list(zip(lookup["chrom"], lookup["strand"], utr0))
    return pd.Series([best.get(k, np.nan) for k in keys], index=lookup.index)


def window_wl_by_transcript(run_dir, id_col):
    """{transcript: window start} of a run made with per-transcript windows (per_transcript.csv
    window_wl, e.g. --cut-at-exon-start), looked up in run_dir or in the run it was derived from
    (params.json source_dir); {} for runs with one window for all transcripts."""
    for d in (run_dir, (read_run_params(run_dir) or {}).get("source_dir")):
        path = os.path.join(d, "per_transcript.csv") if d else None
        if path and os.path.exists(path):
            pt = pd.read_csv(path)
            if "window_wl" in pt.columns:
                return dict(zip(pt[id_col], pt["window_wl"].astype(int)))
    return {}


def anchor_of_run(run_dir):
    """Anchor a results directory was produced with (runs predating --anchor used the PAS)."""
    return (read_run_params(run_dir) or {}).get("anchor", "pas")


# ---------------------------------------------------------------------------
# m6A sites, indexed per chrom+strand for fast window lookup
# ---------------------------------------------------------------------------
def load_m6a_sites(bed_path):
    cols = ["chrom", "start", "end", "site_id", "score", "strand"] + [f"c{i}" for i in range(23)]
    m6a = pd.read_csv(bed_path, sep="\t", header=None, names=cols, usecols=range(6))
    by_key = {}
    for (chrom, strand), g in m6a.groupby(["chrom", "strand"]):
        order = np.argsort(g["start"].values)
        by_key[(chrom, strand)] = (g["start"].values[order], g["site_id"].values[order])
    return by_key


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def window_genomic_bounds(strand, pas, wl=WL, wr=WR):
    """0-based half-open [start, end) forward-strand genomic bounds of the [wl, wr) window."""
    if strand == "+":
        return pas + wl, pas + wr
    return pas - wr, pas - wl


def fetch_window(fasta, chrom, strand, pas, wl=WL, wr=WR):
    """Return the [wl:wr) nt window around the PAS, in 5'->3' transcript orientation,
    plus the 0-based genomic start used for the fetch on the *forward* strand
    (needed later to map m6A genomic coordinates back into window indices)."""
    g_start, g_end = window_genomic_bounds(strand, pas, wl, wr)
    if chrom not in fasta.references:
        return None, None
    if g_start < 0 or g_end > fasta.get_reference_length(chrom):
        return None, None
    seq = fasta.fetch(chrom, g_start, g_end).upper()
    if len(seq) != wr - wl:
        return None, None
    if strand == "-":
        seq = revcomp(seq)
    seq = seq.replace("T", "U")
    return seq, g_start


def genomic_to_window_index(g0, strand, pas, wl=WL):
    """0-based genomic position -> 1-based RNAfold index within the extracted window."""
    if strand == "+":
        i = g0 - (pas + wl) + 1
    else:
        i = pas - wl - g0
    return i


def window_index_to_genomic(i, strand, pas, wl=WL):
    """Inverse of genomic_to_window_index: 1-based RNAfold index -> 0-based genomic
    position. Used to map the switch-competent partner j (and its U-tract flank)
    back to genomic coordinates so it can be intersected with CLIP crosslink sites."""
    if strand == "+":
        g0 = pas + wl + i - 1
    else:
        g0 = pas - wl - i
    return g0


# ---------------------------------------------------------------------------
# Folding / switch-competence
# ---------------------------------------------------------------------------
def fold_bpp(seq, max_bp_span=None):
    """Fold `seq` and return, per 1-based position, the most probable pairing
    partner and its probability. Builds the symmetric bpp matrix via one
    vectorized slice-assignment per row (not a Python-level O(n^2) double
    loop) since this runs once per transcript, thousands of times.
    `max_bp_span` (nt) caps the distance between paired bases; None = unrestricted
    (the setting of every PAS-anchored run)."""
    if max_bp_span is None:
        fc = RNA.fold_compound(seq)
    else:
        md = RNA.md()
        md.max_bp_span = max_bp_span
        fc = RNA.fold_compound(seq, md)
    _, mfe_en = fc.mfe()
    fc.exp_params_rescale(mfe_en)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    mat = np.zeros((n + 1, n + 1), dtype=float)
    for i in range(1, n + 1):
        row = bpp[i]
        L = len(row)
        if L > i:
            mat[i, i:L] = row[i:L]
    mat += mat.T
    best_j = np.argmax(mat, axis=1)
    best_p = mat[np.arange(n + 1), best_j]
    return best_j, best_p


def u_tract_bounds(seq_1indexed, j, n):
    """1-based inclusive (left, right) bounds of the maximal run of consecutive U's that
    contains position j, extending left/right from j itself (no fixed flank). Returns
    (None, None) if position j is not a U."""
    if seq_1indexed[j] != "U":
        return None, None
    left = j
    while left > 1 and seq_1indexed[left - 1] == "U":
        left -= 1
    right = j
    while right < n and seq_1indexed[right + 1] == "U":
        right += 1
    return left, right


def has_u_tract(seq_1indexed, j, n, min_run=U_TRACT_MIN):
    """Is position j itself part of a run of >= min_run consecutive U's?"""
    left, right = u_tract_bounds(seq_1indexed, j, n)
    if left is None:
        return False
    return (right - left + 1) >= min_run


# ---------------------------------------------------------------------------
# Isoform selection (shared by every script in this pipeline)
# ---------------------------------------------------------------------------
ISOFORMS = ("proximal", "distal")


def id_columns(isoform):
    """(key column, partner column) for an isoform: proximal -> (proximal_id, distal_id)."""
    if isoform not in ISOFORMS:
        raise ValueError(f"isoform must be one of {ISOFORMS}, got {isoform!r}")
    return ("proximal_id", "distal_id") if isoform == "proximal" else ("distal_id", "proximal_id")


def load_isoform_table(isoform):
    """One row per isoform of the requested type, with its gene, its partner isoform in the
    pair, and the pair's nuclear-lengthened label. The pair set is table_of_pairs.csv reduced
    to the longest distal per proximal; the label is "distal_id is in SHIFT_FILE". For distal,
    that same pair set is keyed by distal_id instead (identical labels pair by pair). Row order
    for proximal is exactly the original script's, so legacy outputs reproduce byte-for-byte."""
    id_col, partner_col = id_columns(isoform)
    table_of_pairs = pd.read_csv(TABLE_OF_PAIRS, usecols=["proximal_id", "gene", "distal_id", "length"])
    pairs = (
        table_of_pairs.sort_values("length", ascending=False)
        .drop_duplicates(subset="proximal_id", keep="first")
        .reset_index(drop=True)
    )
    shift_df = pd.read_csv(SHIFT_FILE)
    shift_df = shift_df[shift_df["distal_id"].isin(pairs["distal_id"])]
    pairs["nuclear_lengthened"] = pairs["distal_id"].isin(set(shift_df["distal_id"]))

    table = pairs.drop_duplicates(subset=id_col, keep="first").reset_index(drop=True)
    print(f"  {len(table_of_pairs)} pairs -> {len(table)} unique {isoform} isoforms "
          f"({table['nuclear_lengthened'].sum()} nuclear-lengthened / {(~table['nuclear_lengthened']).sum()} other)")
    return table[[id_col, "gene", partner_col, "nuclear_lengthened"]]


# ---------------------------------------------------------------------------
# Provenance: every run writes params.json next to its outputs
# ---------------------------------------------------------------------------
_PROVENANCE_PACKAGES = ["numpy", "pandas", "scipy", "pysam", "ViennaRNA"]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_params(out_dir, run_params, script_paths, input_paths, filename="params.json"):
    """Record what produced a results directory: run parameters, input file paths,
    sha256 of the scripts that ran (the git commit alone doesn't pin them -- these scripts
    may be uncommitted), git commit/dirty state, package versions, and the command line."""
    repo_dir = str(REPO_ROOT)

    def git(*args):
        try:
            return subprocess.run(["git", "-C", repo_dir, *args], capture_output=True, text=True,
                                  check=False).stdout.strip()
        except OSError:
            return None

    versions = {}
    for pkg in _PROVENANCE_PACKAGES:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = None

    params = {
        "run": run_params,
        "inputs": input_paths,
        "scripts_sha256": {os.path.basename(p): _sha256(p) for p in script_paths},
        "git_commit": git("rev-parse", "HEAD"),
        "git_scripts_dirty": bool(git("status", "--porcelain", "--", "src/utrals/m6a_switch", "scripts/m6a_switch")),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": versions,
        "command": [sys.executable, *sys.argv],
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    with open(os.path.join(out_dir, filename), "w") as f:
        json.dump(params, f, indent=2)


def read_run_params(run_dir):
    """The "run" block of a results directory's params.json, or None for legacy directories
    produced before params.json existed."""
    path = os.path.join(run_dir, "params.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)["run"]


def resolve_window(source_dir, isoform, wl, wr):
    """isoform/wl/wr for a downstream step: taken from the upstream run's params.json when it
    exists (so a step can't silently run on a different window than its input), else from the
    command line (legacy directories). Raises if both are given and disagree."""
    upstream = read_run_params(source_dir)
    if upstream is None:
        return isoform, wl, wr
    up = (upstream["isoform"], upstream["wl"], upstream["wr"])
    given = (isoform, wl, wr)
    if any(g is not None and g != u for g, u in zip(given, up)):
        raise ValueError(f"--isoform/--wl/--wr {given} disagree with {source_dir}/params.json {up}")
    return up


# ---------------------------------------------------------------------------
# Per-transcript worker (module-level so multiprocessing can use it; state is inherited from
# the parent via fork, and each worker opens its own FastaFile handle)
# ---------------------------------------------------------------------------
_WORKER = {}


def _init_worker():
    _WORKER["fasta"] = pysam.FastaFile(GENOME_FASTA)


def _process_transcript(task):
    tx_id, gene, partner, nuclear_lengthened, chrom, strand, pas, wl = task   # wl: per-transcript window start
    wr, id_col, partner_col = _WORKER["wr"], _WORKER["id_col"], _WORKER["partner_col"]

    seq, _ = fetch_window(_WORKER["fasta"], chrom, strand, pas, wl, wr)
    if seq is None:
        return None

    per_site_rows = []
    n_sites_in_window = 0
    n_switch_competent = 0

    sites = _WORKER["m6a_by_key"].get((chrom, strand))
    if sites is not None:
        starts, ids = sites
        g_lo, g_hi = window_genomic_bounds(strand, pas, wl, wr)
        lo_idx = np.searchsorted(starts, g_lo, side="left")
        hi_idx = np.searchsorted(starts, g_hi, side="left")
        if hi_idx > lo_idx:
            best_j, best_p = fold_bpp(seq)
            seq_1idx = [""] + list(seq)  # 1-indexed access
            n = len(seq)
            for g0, site_id in zip(starts[lo_idx:hi_idx], ids[lo_idx:hi_idx]):
                i = genomic_to_window_index(g0, strand, pas, wl)
                if i < 1 or i > n:
                    continue
                j = int(best_j[i])
                p = float(best_p[i])
                u_tract = has_u_tract(seq_1idx, j, n) if j > 0 else False
                switch_competent = bool(j > 0 and p >= P_THRESHOLD and u_tract)
                n_sites_in_window += 1
                n_switch_competent += int(switch_competent)

                if j > 0:
                    g_j = window_index_to_genomic(j, strand, pas, wl)
                    left, right = u_tract_bounds(seq_1idx, j, n)
                    if left is not None:
                        g_left = window_index_to_genomic(left, strand, pas, wl)
                        g_right = window_index_to_genomic(right, strand, pas, wl)
                        # genomic coordinates increase with window index on '+' strand and
                        # decrease on '-' strand -- normalize to a sorted (start, end) pair
                        u_tract_g_start, u_tract_g_end = sorted((g_left, g_right))
                    else:
                        u_tract_g_start = u_tract_g_end = np.nan
                else:
                    g_j = u_tract_g_start = u_tract_g_end = np.nan

                per_site_rows.append({
                    id_col: tx_id,
                    "gene": gene,
                    "chrom": chrom,
                    "strand": strand,
                    "site_id": site_id,
                    "i": i,
                    "j": j,
                    "p_ij": p,
                    "u_tract": u_tract,
                    "switch_competent": switch_competent,
                    "j_genomic": g_j,
                    "u_tract_genomic_start": u_tract_g_start,
                    "u_tract_genomic_end": u_tract_g_end,
                    "nuclear_lengthened": nuclear_lengthened,
                })

    per_tx_row = {
        id_col: tx_id,
        "gene": gene,
        partner_col: partner,
        "chrom": chrom,
        "strand": strand,
        "n_m6a_sites": n_sites_in_window,
        "n_switch_competent": n_switch_competent,
        "frac_switch_competent": (n_switch_competent / n_sites_in_window) if n_sites_in_window > 0 else np.nan,
        "nuclear_lengthened": nuclear_lengthened,
        "n_A": seq.count("A"),
        "window_wl": wl,
    }
    return per_site_rows, per_tx_row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--isoform", choices=ISOFORMS, default="proximal",
                    help="which isoform of each pair to analyse (default: proximal)")
    ap.add_argument("--anchor", choices=ANCHORS, default="pas",
                    help="window anchor: the PAS (default) or the 3'UTR start (first nt after the stop codon)")
    ap.add_argument("--cut-at-exon-start", action="store_true",
                    help="with --anchor utr_start: start the window at the first nt of the exon holding the "
                         "stop codon when that exon begins inside [wl:0) (no intron folded)")
    ap.add_argument("--wl", type=int, default=WL,
                    help=f"window start relative to the PAS, transcript orientation (default: {WL})")
    ap.add_argument("--wr", type=int, default=WR,
                    help=f"window end (exclusive) relative to the PAS (default: {WR})")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: legacy m6A_switch_competence_<wl>_<wr> for "
                         "proximal; required for distal)")
    ap.add_argument("--n-jobs", type=int, default=1,
                    help="worker processes for folding (output is identical for any value)")
    args = ap.parse_args(argv)
    if args.wr <= args.wl:
        ap.error("--wr must be greater than --wl")
    if args.cut_at_exon_start and args.anchor != "utr_start":
        ap.error("--cut-at-exon-start only applies to --anchor utr_start")
    if args.out_dir is None:
        if args.isoform != "proximal":
            ap.error("--out-dir is required for --isoform distal")
        args.out_dir = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_{args.wl}_{args.wr}"
    return args


def main(argv=None):
    args = parse_args(argv)
    id_col, partner_col = id_columns(args.isoform)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"isoform={args.isoform} anchor={args.anchor} window=[{args.wl}:{args.wr}) out_dir={args.out_dir} n_jobs={args.n_jobs}")

    print("Loading table_of_pairs.csv + shift file ...")
    tx_table = load_isoform_table(args.isoform)

    print(f"Parsing GTF for PAS coordinates (pre-filtered to the {args.isoform} isoforms we need) ...")
    needed_ids = set(tx_table[id_col])
    pas_lookup = load_pas_lookup(GTF_PATH, needed_ids=needed_ids)
    print(f"  {len(pas_lookup)} / {len(needed_ids)} needed transcripts found with coordinates")

    print("Loading m6A sites ...")
    _WORKER.update(m6a_by_key=load_m6a_sites(M6A_BED), wr=args.wr,
                   id_col=id_col, partner_col=partner_col)

    wl_by_tx = pd.Series(args.wl, index=pas_lookup.index)
    n_cut = 0
    if args.cut_at_exon_start:
        print("Matching stop codons to GENCODE exons to cut windows at the stop-codon exon start ...")
        nt_before = stop_exon_nt_before_utr(pas_lookup)
        cut = nt_before.notna() & (-nt_before > args.wl)
        wl_by_tx[cut] = -nt_before[cut].astype(int)
        n_cut = int(cut.sum())
        print(f"  {nt_before.notna().sum()} / {len(nt_before)} matched; {n_cut} windows cut at the exon start")

    tasks = []
    n_no_coords = 0
    for row in tx_table.itertuples(index=False):
        tx_id = getattr(row, id_col)
        if tx_id not in pas_lookup.index:
            n_no_coords += 1
            continue
        info = pas_lookup.loc[tx_id]
        tasks.append((tx_id, row.gene, getattr(row, partner_col), row.nuclear_lengthened,
                      info["chrom"], info["strand"], anchor_coordinate(info, args.anchor), int(wl_by_tx[tx_id])))

    if args.n_jobs > 1:
        pool = Pool(args.n_jobs, initializer=_init_worker)
        results = pool.imap(_process_transcript, tasks, chunksize=8)  # imap preserves input order
    else:
        pool = None
        _init_worker()
        results = map(_process_transcript, tasks)

    per_site_rows = []
    per_tx_rows = []
    n_no_seq = 0
    n_ok = 0
    for result in results:
        if result is None:
            n_no_seq += 1
            continue
        n_ok += 1
        site_rows, tx_row = result
        per_site_rows.extend(site_rows)
        per_tx_rows.append(tx_row)
        if n_ok % 1000 == 0:
            print(f"  processed {n_ok} transcripts ...", flush=True)
    if pool is not None:
        pool.close()
        pool.join()

    print(f"Done. ok={n_ok}, no_coords={n_no_coords}, no_seq={n_no_seq}")

    per_site_df = pd.DataFrame(per_site_rows)
    per_tx_df = pd.DataFrame(per_tx_rows)

    per_site_path = os.path.join(args.out_dir, "per_site.csv")
    per_tx_path = os.path.join(args.out_dir, "per_transcript.csv")
    per_site_df.to_csv(per_site_path, index=False)
    per_tx_df.to_csv(per_tx_path, index=False)
    print(f"Wrote {per_site_path} ({len(per_site_df)} rows)")
    print(f"Wrote {per_tx_path} ({len(per_tx_df)} rows)")

    write_params(
        args.out_dir,
        run_params={"step": "switch_competence_u_tract", "isoform": args.isoform, "anchor": args.anchor,
                    "cut_at_exon_start": args.cut_at_exon_start, "n_windows_cut": n_cut, "wl": args.wl,
                    "wr": args.wr, "p_threshold": P_THRESHOLD, "u_tract_min": U_TRACT_MIN,
                    "n_transcripts": n_ok, "n_no_coords": n_no_coords, "n_no_seq": n_no_seq},
        script_paths=[os.path.abspath(__file__)],
        input_paths={"table_of_pairs": TABLE_OF_PAIRS, "shift_file": SHIFT_FILE, "gtf": GTF_PATH,
                     "m6a_bed": M6A_BED, "genome_fasta": GENOME_FASTA,
                     **({"gencode_gtf_exons": GENCODE_GTF} if args.cut_at_exon_start else {})},
    )


if __name__ == "__main__":
    main()
