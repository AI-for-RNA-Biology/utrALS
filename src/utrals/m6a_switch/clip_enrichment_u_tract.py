#!/usr/bin/env python
"""
Do the candidate RBPs in --candidates bind specifically at the U-tracts of switch-competent
m6A sites (identified by scripts/m6a_switch/m6a_switch_competence_analysis.py), or do they just bind
U-tracts in general (many of these RBPs -- HuR, TIAL1, HNRNPC, CELF4 -- are known
U-rich-element binders regardless of any m6A structural-switching context)?

To isolate the m6A-switching effect from plain U-richness preference, the background here is
NOT "all nucleotides in the window" -- it's specifically *other* U-tracts in the same window
that are not linked to a switch-competent m6A site. Both the test set and the background are
real, maximal U-runs (>=3 consecutive U, scanned directly from the genome), so any RBP's
baseline affinity for U-tracts as such is matched between the two groups; only their
association with a switch-competent m6A site differs.

The isoform (proximal/distal) and window are those of the switch-competence run given by
--switch-dir (read from its params.json; --isoform/--wl/--wr only needed for legacy
directories that predate params.json). The scan window is always identical to the fold
window. Run with no arguments to reproduce the original proximal [0:200] run.

Note: the R2 regulators of lengthening (default --candidates) were selected on the proximal [0,100]
"region2" window (Figure_regulators_lengthening.ipynb); applying it to other windows/isoforms asks whether those
same RBPs show the preference there.

This does NOT re-run RNAfold -- it reuses the switch-competence run's per_site.csv directly,
including its u_tract_genomic_start/end columns (0-based, inclusive: the actual
consecutive-U run containing j).

Pipeline:
  1. For every isoform, scan the window for all maximal U-runs (>=3 nt).
  2. Label each U-run "switch" if it overlaps the U-tract run of a switch-competent m6A site
     in that transcript, else "other".
  3. For each RBP: n = crosslinks landing in ANY U-run (switch or other), k = of those, how
     many land in a "switch" U-run, p0 = (switch U-run nt) / (all U-run nt). One-sided
     binomial test: is k/n > p0?

Outputs under --out-dir:
  - u_tract_catalog.csv : one row per U-run in the window (all transcripts), switch vs. other
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
LEGACY_SWITCH_DIR = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_0_200"
LEGACY_OUT_DIR = f"{M6A_RESULTS_ROOT}/clip_u_tract_enrichment_0_200"
CANDIDATES_CSV = (f"{REPO_ROOT}/Tables/Figure_regulators_lengthening/"
                  "regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_region2_downstream_pPAS_APA_vs_nuclear_retention_class.csv")
CLIP_ROOT = CLIP_DIR

U_RUN_MIN = 3


# ---------------------------------------------------------------------------
# scan-window genomic coordinate mapping (delegates to the switch-competence module so the
# scan window and the fold window can never disagree)
# ---------------------------------------------------------------------------
def local_run_to_genomic(local_start, local_end, strand, pas, wl):
    """0-based half-open local run [local_start, local_end) -> genomic half-open interval."""
    g_first = window_index_to_genomic(local_start + 1, strand, pas, wl)
    g_last = window_index_to_genomic(local_end, strand, pas, wl)
    g_lo, g_hi = sorted((g_first, g_last))
    return g_lo, g_hi + 1


def find_u_runs(seq, min_run=U_RUN_MIN):
    """Maximal runs of >=min_run consecutive U's. Returns [(local_start, local_end), ...]."""
    runs = []
    i, n = 0, len(seq)
    while i < n:
        if seq[i] == "U":
            j = i
            while j < n and seq[j] == "U":
                j += 1
            if j - i >= min_run:
                runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


# ---------------------------------------------------------------------------
# Switch-competent U-tract spans (from the switch-competence run's per_site.csv --
# u_tract_genomic_start/end are 0-based INCLUSIVE, converted here to half-open)
# ---------------------------------------------------------------------------
def switch_flanks_by_transcript(switch_dir, id_col):
    per_site = pd.read_csv(os.path.join(switch_dir, "per_site.csv"))
    switch = per_site[per_site["switch_competent"]].copy()
    print(f"{len(switch)} switch-competent sites in {switch_dir}")

    flanks = {}
    for r in switch.itertuples(index=False):
        f_start, f_end = int(r.u_tract_genomic_start), int(r.u_tract_genomic_end) + 1
        flanks.setdefault(getattr(r, id_col), []).append((f_start, f_end))
    return flanks


# ---------------------------------------------------------------------------
# Build the full U-tract catalog (all isoforms), labeled switch vs. other
# ---------------------------------------------------------------------------
def build_u_tract_catalog(switch_dir, isoform, wl, wr):
    id_col, _ = id_columns(isoform)
    flanks_by_tx = switch_flanks_by_transcript(switch_dir, id_col)
    wl_of = window_wl_by_transcript(switch_dir, id_col)   # per-transcript window starts, if any

    tx_table = load_isoform_table(isoform)
    pas_lookup = load_pas_lookup(GTF_PATH, needed_ids=set(tx_table[id_col]))
    fasta = pysam.FastaFile(GENOME_FASTA)

    rows = []
    n_done = 0
    for r in tx_table.itertuples(index=False):
        tx_id = getattr(r, id_col)
        n_done += 1
        if n_done % 2000 == 0:
            print(f"  scanned {n_done} / {len(tx_table)} transcripts for U-runs ...")
        if tx_id not in pas_lookup.index:
            continue
        info = pas_lookup.loc[tx_id]
        chrom, strand, pas = info["chrom"], info["strand"], anchor_coordinate(info, anchor_of_run(switch_dir))

        wl_t = wl_of.get(tx_id, wl)
        seq, _ = fetch_window(fasta, chrom, strand, pas, wl_t, wr)
        if seq is None:
            continue

        runs = find_u_runs(seq)
        if not runs:
            continue

        flanks = flanks_by_tx.get(tx_id, [])
        for local_start, local_end in runs:
            g_start, g_end = local_run_to_genomic(local_start, local_end, strand, pas, wl_t)
            is_switch = any(g_start < f_end and f_start < g_end for f_start, f_end in flanks)
            rows.append({
                id_col: tx_id, "gene": r.gene, "chrom": chrom, "strand": strand,
                "u_run_genomic_start": g_start, "u_run_genomic_end": g_end,
                "u_run_len": g_end - g_start, "category": "switch" if is_switch else "other",
            })

    df = pd.DataFrame(rows)
    print(f"{len(df)} U-runs (>= {U_RUN_MIN} nt) found in [{wl}:{wr}] across "
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


def test_rbp(clip_by_key, tx_runs):
    """tx_runs: {(chrom, strand): [(start, end, category), ...]} grouped for fast lookup
    via a parallel sorted-starts array per (chrom, strand)."""
    n_total = 0
    k_switch = 0
    for (chrom, strand), (starts, ends, categories) in tx_runs.items():
        sites = clip_by_key.get((chrom, strand))
        if sites is None:
            continue
        # for each U-run, count crosslinks landing inside it
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
                    help="U-tract switch-competence run (default: legacy proximal [0:200])")
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

    catalog = build_u_tract_catalog(args.switch_dir, args.isoform, args.wl, args.wr)
    catalog.to_csv(os.path.join(args.out_dir, "u_tract_catalog.csv"), index=False)

    total_switch_nt = catalog.loc[catalog["category"] == "switch", "u_run_len"].sum()
    total_other_nt = catalog.loc[catalog["category"] == "other", "u_run_len"].sum()
    p0 = total_switch_nt / (total_switch_nt + total_other_nt)
    print(f"Background: {total_switch_nt} switch-associated U-run nt / "
          f"{total_switch_nt + total_other_nt} total U-run nt in [{args.wl}:{args.wr}] -> p0={p0:.4f}")

    # group runs by (chrom, strand) for fast per-RBP lookup
    tx_runs = {}
    for (chrom, strand), g in catalog.groupby(["chrom", "strand"]):
        g = g.sort_values("u_run_genomic_start")
        tx_runs[(chrom, strand)] = (
            g["u_run_genomic_start"].values, g["u_run_genomic_end"].values, g["category"].values,
        )

    if args.catalog_only:
        _write_params(args, p0)
        print(f"catalog only: wrote {os.path.join(args.out_dir, 'u_tract_catalog.csv')}")
        return

    candidates = pd.read_csv(args.candidates, index_col=0)
    if "CLIP" not in candidates.columns and "regulator" in candidates.columns:
        candidates = candidates.rename(columns={"regulator": "CLIP"})  # newer candidate tables

    results = []
    for r in candidates.itertuples(index=False):
        print(f"Processing {r.RBP} ({r.CLIP}) ...")
        gff_path = resolve_clip_file(r.CLIP)
        clip_by_key = load_crosslinks(gff_path)
        n, k = test_rbp(clip_by_key, tx_runs)
        if n == 0:
            results.append({"CLIP": r.CLIP, "RBP": r.RBP, "category": r.category,
                             "n_crosslinks_in_u_tracts": 0, "n_in_switch_u_tract": 0,
                             "observed_rate": np.nan, "expected_rate": p0,
                             "odds_ratio": np.nan, "p_value": np.nan})
            continue
        test = binomtest(k, n, p0, alternative="greater")
        obs_rate = k / n
        odds_ratio = (obs_rate / (1 - obs_rate)) / (p0 / (1 - p0)) if 0 < obs_rate < 1 else np.nan
        results.append({
            "CLIP": r.CLIP, "RBP": r.RBP, "category": r.category,
            "n_crosslinks_in_u_tracts": n, "n_in_switch_u_tract": k,
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
        run_params={"step": "clip_u_tract_enrichment", "isoform": args.isoform, "anchor": anchor_of_run(args.switch_dir), "wl": args.wl,
                    "wr": args.wr, "u_run_min": U_RUN_MIN, "switch_dir": os.path.abspath(args.switch_dir),
                    "p0": float(p0), "catalog_only": args.catalog_only},
        script_paths=[os.path.abspath(__file__),
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), "competence.py")],
        input_paths={"candidates": None if args.catalog_only else os.path.abspath(args.candidates), "clip_root": CLIP_ROOT,
                     "gtf": GTF_PATH, "genome_fasta": GENOME_FASTA},
    )


if __name__ == "__main__":
    main()
