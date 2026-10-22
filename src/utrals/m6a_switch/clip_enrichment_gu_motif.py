#!/usr/bin/env python
"""
GU-motif analogue of clip_u_tract_switch_competent_enrichment.py: do the candidate RBPs in
--candidates bind preferentially at UGU/UUG/GUG/GUU 3-mers that are part of a switch-competent
m6A site (per scripts/m6a_switch/m6a_switch_competence_guu_motif_analysis.py), or do they just bind these
motifs in general?

Same logic as the U-tract version: the background is NOT "all nucleotides" -- it's other
occurrences of the same four motifs that are *not* linked to a switch-competent m6A site.
Both the test set and the background are real, exact 3-mer matches scanned directly from the
genome, so an RBP's baseline affinity for GU-rich 3-mers as such is matched between them;
only their association with a switch-competent m6A site differs.

One difference from the U-tract version: individual 3-mer occurrences of UGU/UUG/GUG/GUU
overlap each other by construction (e.g. "GUGUU" contains GUG, UGU, and GUU as overlapping
matches), unlike maximal U-runs which by definition don't overlap. Rather than merging
overlapping occurrences into runs (which would blur which specific motif identity is
switch-linked), this catalog keeps every individual 3-mer occurrence as its own entry, and
uses occurrence COUNTS (not nucleotide spans) as the background unit, since every entry has
the same fixed length. This means a single genomic position covered by two overlapping
occurrences of different categories can contribute to both "switch" and "other" counts, and
a single CLIP crosslink landing in that overlap is counted once per occurrence it falls in --
a modeling simplification, not a bug; flagged here rather than silently glossed over.

The isoform (proximal/distal) and window are those of the GU-motif switch-competence run given
by --switch-dir (read from its params.json; --isoform/--wl/--wr only needed for legacy
directories). Run with no arguments to reproduce the original proximal [0:200] run.

Pipeline:
  1. For every isoform, scan the window (identical to the fold window) for every exact
     occurrence of UGU/UUG/GUG/GUU.
  2. Label each occurrence "switch" if it IS one of the specific 3-mer window(s) that made a
     real m6A site switch-competent in that transcript, else "other".
  3. For each RBP: n = crosslinks landing in ANY motif occurrence (switch or other), k = of
     those, how many land in a "switch" occurrence, p0 = (switch occurrences) / (all
     occurrences). One-sided binomial test: is k/n > p0?

Outputs under --out-dir:
  - motif_catalog.csv   : one row per motif occurrence in the window (all transcripts),
                           switch vs. other
  - clip_enrichment.csv : one row per RBP, with n, k, p0, observed rate, OR, p-value
  - params.json
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import pysam
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from ..paths import CLIP_DIR, REPO_ROOT  # data locations, see paths.env
from .competence import (
    GTF_PATH, GENOME_FASTA, M6A_RESULTS_ROOT, ISOFORMS, load_pas_lookup, fetch_window, anchor_coordinate, anchor_of_run, window_wl_by_transcript,
    window_index_to_genomic, id_columns, load_isoform_table, resolve_window, write_params,
)

# ---------------------------------------------------------------------------
# Defaults (legacy proximal [0:200] run)
# ---------------------------------------------------------------------------
LEGACY_SWITCH_DIR = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_0_200_guu_motif"
LEGACY_OUT_DIR = f"{M6A_RESULTS_ROOT}/clip_guu_motif_enrichment_0_200"
CANDIDATES_CSV = (f"{REPO_ROOT}/Tables/Figure_regulators_lengthening/"
                  "regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_region2_downstream_pPAS_APA_vs_nuclear_retention_class.csv")
CLIP_ROOT = CLIP_DIR

MOTIFS = ["UGU", "UUG", "GUG", "GUU"]
MOTIF_LEN = 3


# ---------------------------------------------------------------------------
# scan-window genomic coordinate mapping (delegates to the switch-competence module so the
# scan window and the fold window can never disagree)
# ---------------------------------------------------------------------------
def local_span_to_genomic(local_start, local_end, strand, pas, wl):
    """0-based half-open local span [local_start, local_end) -> genomic half-open interval."""
    g_first = window_index_to_genomic(local_start + 1, strand, pas, wl)
    g_last = window_index_to_genomic(local_end, strand, pas, wl)
    g_lo, g_hi = sorted((g_first, g_last))
    return g_lo, g_hi + 1


def find_motif_occurrences(seq):
    """Every exact occurrence of any MOTIFS string, sliding one nt at a time (occurrences
    overlap by construction -- not merged). Returns [(local_start, local_end, motif), ...]."""
    occurrences = []
    n = len(seq)
    for start in range(n - MOTIF_LEN + 1):
        sub = seq[start:start + MOTIF_LEN]
        if sub in MOTIFS:
            occurrences.append((start, start + MOTIF_LEN, sub))
    return occurrences


# ---------------------------------------------------------------------------
# Switch-competent motif occurrences (from the GU-motif run's per_site.csv; recompute the
# exact matching 3-mer window(s) at j from the real sequence -- per_site.csv only stores
# boolean motif flags, not their genomic position)
# ---------------------------------------------------------------------------
def switch_occurrences_by_transcript(switch_dir, id_col, wl, wr):
    per_site = pd.read_csv(os.path.join(switch_dir, "per_site.csv"))
    switch = per_site[per_site["switch_competent"]].copy()
    print(f"{len(switch)} switch-competent sites in {switch_dir}")

    pas_lookup = load_pas_lookup(GTF_PATH, needed_ids=set(switch[id_col]))
    wl_of = window_wl_by_transcript(switch_dir, id_col)   # per-transcript window starts, if any
    fasta = pysam.FastaFile(GENOME_FASTA)

    seq_cache = {}
    occurrences = {}
    for r in switch.itertuples(index=False):
        tx_id = getattr(r, id_col)
        info = pas_lookup.loc[tx_id]
        strand, pas = info["strand"], anchor_coordinate(info, anchor_of_run(switch_dir))
        if tx_id not in seq_cache:
            seq_cache[tx_id], _ = fetch_window(fasta, info["chrom"], strand, pas, wl_of.get(tx_id, wl), wr)
        seq = seq_cache[tx_id]
        if seq is None:
            continue

        j = int(r.j)  # 1-based
        n = len(seq)
        for offset in range(MOTIF_LEN):
            start0 = (j - 1) - offset  # 0-based local start of the 3-mer covering j
            end0 = start0 + MOTIF_LEN
            if start0 < 0 or end0 > n:
                continue
            if seq[start0:end0] not in MOTIFS:
                continue
            occurrences.setdefault(tx_id, []).append(local_span_to_genomic(start0, end0, strand, pas, wl_of.get(tx_id, wl)))
    return occurrences


# ---------------------------------------------------------------------------
# Build the full motif catalog (all isoforms), labeled switch vs. other
# ---------------------------------------------------------------------------
def build_motif_catalog(switch_dir, isoform, wl, wr):
    id_col, _ = id_columns(isoform)
    switch_by_tx = switch_occurrences_by_transcript(switch_dir, id_col, wl, wr)

    tx_table = load_isoform_table(isoform)
    pas_lookup = load_pas_lookup(GTF_PATH, needed_ids=set(tx_table[id_col]))
    wl_of = window_wl_by_transcript(switch_dir, id_col)   # per-transcript window starts, if any
    fasta = pysam.FastaFile(GENOME_FASTA)

    rows = []
    n_done = 0
    for r in tx_table.itertuples(index=False):
        tx_id = getattr(r, id_col)
        n_done += 1
        if n_done % 2000 == 0:
            print(f"  scanned {n_done} / {len(tx_table)} transcripts for motif occurrences ...")
        if tx_id not in pas_lookup.index:
            continue
        info = pas_lookup.loc[tx_id]
        chrom, strand, pas = info["chrom"], info["strand"], anchor_coordinate(info, anchor_of_run(switch_dir))

        wl_t = wl_of.get(tx_id, wl)
        seq, _ = fetch_window(fasta, chrom, strand, pas, wl_t, wr)
        if seq is None:
            continue

        occurrences = find_motif_occurrences(seq)
        if not occurrences:
            continue

        switch_spans = switch_by_tx.get(tx_id, [])
        for local_start, local_end, motif in occurrences:
            g_start, g_end = local_span_to_genomic(local_start, local_end, strand, pas, wl_t)
            is_switch = any(g_start < s_end and s_start < g_end for s_start, s_end in switch_spans)
            rows.append({
                id_col: tx_id, "gene": r.gene, "chrom": chrom, "strand": strand,
                "motif": motif, "occurrence_genomic_start": g_start, "occurrence_genomic_end": g_end,
                "category": "switch" if is_switch else "other",
            })

    df = pd.DataFrame(rows)
    print(f"{len(df)} motif occurrences found in [{wl}:{wr}] across "
          f"{df[id_col].nunique()} transcripts "
          f"({(df['category'] == 'switch').sum()} switch-associated, "
          f"{(df['category'] == 'other').sum()} other)")
    return df


# ---------------------------------------------------------------------------
# CLIP crosslinks
# ---------------------------------------------------------------------------
def resolve_clip_file(clip_id):
    """Candidate table 'CLIP' identifier -> crosslink gff path."""
    if clip_id.endswith("_iCLIP"):
        substr = clip_id[: -len("_iCLIP")]
        matches = [f for f in glob.glob(os.path.join(CLIP_ROOT, "iCLIP", "*.gff"))
                   if substr.lower() in os.path.basename(f).lower() and not f.endswith(".unmap")]
    else:
        prefix, source = clip_id.rsplit("_", 1)  # e.g. "DDX6_2", "K652"
        matches = [f for f in glob.glob(os.path.join(CLIP_ROOT, "eCLIP", source, f"{prefix}_ENC*.gff"))
                   if not f.endswith(".unmap")]
        non_sorted = [f for f in matches if not f.endswith(".sorted.gff")]
        if non_sorted:
            matches = non_sorted
    if len(matches) != 1:
        raise ValueError(f"Expected exactly 1 file for CLIP id {clip_id!r}, found {len(matches)}: {matches}")
    return matches[0]


def load_crosslinks(gff_path):
    df = pd.read_csv(gff_path, sep="\t", header=None, usecols=[0, 3, 4, 6],
                      names=["chrom", "start", "end", "strand"])
    df["pos0"] = df["start"] - 1  # gff is 1-based inclusive; single-nt features here
    by_key = {}
    for (chrom, strand), g in df.groupby(["chrom", "strand"]):
        by_key[(chrom, strand)] = np.sort(g["pos0"].values)
    return by_key


def test_rbp(clip_by_key, tx_occurrences):
    """tx_occurrences: {(chrom, strand): (starts, ends, categories)} grouped for fast lookup."""
    n_total = 0
    k_switch = 0
    for (chrom, strand), (starts, ends, categories) in tx_occurrences.items():
        sites = clip_by_key.get((chrom, strand))
        if sites is None:
            continue
        lo_idx = np.searchsorted(sites, starts, side="left")
        hi_idx = np.searchsorted(sites, ends, side="left")
        for lo, hi, cat in zip(lo_idx, hi_idx, categories):
            count = hi - lo
            if count <= 0:
                continue
            n_total += count
            if cat == "switch":
                k_switch += count
    return n_total, k_switch


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--switch-dir", default=LEGACY_SWITCH_DIR,
                    help="GU-motif switch-competence run (default: legacy proximal [0:200])")
    ap.add_argument("--out-dir", default=None, help="output directory (default: legacy, proximal only)")
    ap.add_argument("--catalog-only", action="store_true",
                    help="only write the element catalog (switch-linked vs. other) and params.json, no candidate "
                         "test -- e.g. for regions tested only in clip_enrichment_all_datasets.py")
    ap.add_argument("--candidates", default=CANDIDATES_CSV, help="candidate RBP table (CLIP or regulator, RBP, category columns)")
    ap.add_argument("--isoform", choices=ISOFORMS, default=None, help="only for legacy --switch-dir; default proximal")
    ap.add_argument("--wl", type=int, default=None, help="only for legacy --switch-dir; default 0")
    ap.add_argument("--wr", type=int, default=None, help="only for legacy --switch-dir; default 200")
    args = ap.parse_args(argv)

    isoform, wl, wr = resolve_window(args.switch_dir, args.isoform, args.wl, args.wr)
    args.isoform = isoform if isoform is not None else "proximal"
    args.wl = wl if wl is not None else 0
    args.wr = wr if wr is not None else 200
    if args.out_dir is None:
        if (args.isoform, args.wl, args.wr) != ("proximal", 0, 200):
            ap.error("--out-dir is required unless reproducing the legacy proximal [0:200] run")
        args.out_dir = LEGACY_OUT_DIR
    return args


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"isoform={args.isoform} window=[{args.wl}:{args.wr}) switch_dir={args.switch_dir} out_dir={args.out_dir}")

    catalog = build_motif_catalog(args.switch_dir, args.isoform, args.wl, args.wr)
    catalog.to_csv(os.path.join(args.out_dir, "motif_catalog.csv"), index=False)

    n_switch = (catalog["category"] == "switch").sum()
    n_other = (catalog["category"] == "other").sum()
    p0 = n_switch / (n_switch + n_other)
    print(f"Background: {n_switch} switch-associated motif occurrences / "
          f"{n_switch + n_other} total occurrences in [{args.wl}:{args.wr}] -> p0={p0:.4f}")

    tx_occurrences = {}
    for (chrom, strand), g in catalog.groupby(["chrom", "strand"]):
        g = g.sort_values("occurrence_genomic_start")
        tx_occurrences[(chrom, strand)] = (
            g["occurrence_genomic_start"].values, g["occurrence_genomic_end"].values, g["category"].values,
        )

    if args.catalog_only:
        _write_params(args, p0)
        print(f"catalog only: wrote {os.path.join(args.out_dir, 'motif_catalog.csv')}")
        return

    candidates = pd.read_csv(args.candidates, index_col=0)
    if "CLIP" not in candidates.columns and "regulator" in candidates.columns:
        candidates = candidates.rename(columns={"regulator": "CLIP"})  # newer candidate tables

    results = []
    for r in candidates.itertuples(index=False):
        print(f"Processing {r.RBP} ({r.CLIP}) ...")
        gff_path = resolve_clip_file(r.CLIP)
        clip_by_key = load_crosslinks(gff_path)
        n, k = test_rbp(clip_by_key, tx_occurrences)
        if n == 0:
            results.append({"CLIP": r.CLIP, "RBP": r.RBP, "category": r.category,
                             "n_crosslinks_in_motifs": 0, "n_in_switch_motif": 0,
                             "observed_rate": np.nan, "expected_rate": p0,
                             "odds_ratio": np.nan, "p_value": np.nan})
            continue
        test = binomtest(k, n, p0, alternative="greater")
        obs_rate = k / n
        odds_ratio = (obs_rate / (1 - obs_rate)) / (p0 / (1 - p0)) if 0 < obs_rate < 1 else np.nan
        results.append({
            "CLIP": r.CLIP, "RBP": r.RBP, "category": r.category,
            "n_crosslinks_in_motifs": n, "n_in_switch_motif": k,
            "observed_rate": obs_rate, "expected_rate": p0,
            "odds_ratio": odds_ratio, "p_value": test.pvalue,
        })

    results_df = pd.DataFrame(results).sort_values("p_value")
    ok = results_df["p_value"].notna()   # BH over the candidate RBPs of this run
    results_df.loc[ok, "q_value_bh"] = multipletests(results_df.loc[ok, "p_value"], method="fdr_bh")[1]
    out_path = os.path.join(args.out_dir, "clip_enrichment.csv")
    results_df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")
    print(results_df.to_string(index=False))

    _write_params(args, p0)


def _write_params(args, p0):
    write_params(
        args.out_dir,
        run_params={"step": "clip_gu_motif_enrichment", "isoform": args.isoform, "anchor": anchor_of_run(args.switch_dir), "wl": args.wl,
                    "wr": args.wr, "motifs": MOTIFS, "switch_dir": os.path.abspath(args.switch_dir),
                    "p0": float(p0), "catalog_only": args.catalog_only},
        script_paths=[os.path.abspath(__file__),
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), "competence.py")],
        input_paths={"candidates": None if args.catalog_only else os.path.abspath(args.candidates), "clip_root": CLIP_ROOT,
                     "gtf": GTF_PATH, "genome_fasta": GENOME_FASTA},
    )


if __name__ == "__main__":
    main()
