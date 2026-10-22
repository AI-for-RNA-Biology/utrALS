#!/usr/bin/env python
"""
GU-motif variant of a U-tract switch-competence run (scripts/m6a_switch/m6a_switch_competence_analysis.py),
with a different switch-competence criterion.

Reuses that run's per_site.csv directly (--source-dir) -- p(i,j) and the partner j don't
depend on the downstream U-tract/motif criterion, so this does NOT re-run RNAfold. Only the
sequence around j is re-examined, over the same isoform/window as the source run (read from
its params.json; --isoform/--wl/--wr are only needed for legacy source directories that
predate params.json).

Here, a measured m6A site at position i is "switch-competent" if:
  1. p(i, j) >= 0.5 to its most probable partner j (same as the source run), AND
  2. position j is part of a "UGU", "UUG", "GUG", or "GUU" motif overlapping position j
     (any 3-nt window of the sequence containing position j reads exactly one of those four
     strings -- the four permutations of one G and two U's).

Run with no arguments to reproduce the original proximal [0:200] run (legacy defaults).

Outputs under --out-dir:
  - per_site.csv       : one row per (transcript, m6A site), with motif flags/type and
                          switch_competent under this new criterion
  - per_transcript.csv : one row per transcript, with switch-competent counts/fractions,
                          the nuclear-lengthened group label, and n_A if the source run has it
  - params.json
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam

from .competence import (
    GTF_PATH, GENOME_FASTA, M6A_RESULTS_ROOT, P_THRESHOLD, ISOFORMS, load_pas_lookup, fetch_window, anchor_coordinate, anchor_of_run, window_wl_by_transcript,
    id_columns, resolve_window, write_params,
)

LEGACY_SOURCE_DIR = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_0_200"

MOTIFS = ["UGU", "UUG", "GUG", "GUU"]
MOTIF_LEN = 3


def motif_at_position(seq_1idx, j, n):
    """Which (if any) of MOTIFS overlaps position j -- any 3-nt window of the sequence
    containing position j that reads exactly one of the four strings. Returns the set of
    matching motifs (usually 0 or 1; could rarely be >1 for a short window like this)."""
    matches = set()
    for offset in range(MOTIF_LEN):  # 3-mer starting at j-offset covers j when offset in 0..2
        start = j - offset
        end = start + MOTIF_LEN
        if start < 1 or end - 1 > n:
            continue
        sub = "".join(seq_1idx[start:end])
        if sub in MOTIFS:
            matches.add(sub)
    return matches


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", default=LEGACY_SOURCE_DIR,
                    help="U-tract switch-competence run to reuse p(i,j) from (default: legacy proximal [0:200])")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: legacy m6A_switch_competence_<wl>_<wr>_guu_motif for proximal)")
    ap.add_argument("--isoform", choices=ISOFORMS, default=None,
                    help="only needed if --source-dir has no params.json (legacy); default proximal")
    ap.add_argument("--wl", type=int, default=None, help="only needed for legacy --source-dir; default 0")
    ap.add_argument("--wr", type=int, default=None, help="only needed for legacy --source-dir; default 200")
    args = ap.parse_args(argv)

    isoform, wl, wr = resolve_window(args.source_dir, args.isoform, args.wl, args.wr)
    args.isoform = isoform if isoform is not None else "proximal"
    args.wl = wl if wl is not None else 0
    args.wr = wr if wr is not None else 200
    if args.out_dir is None:
        if args.isoform != "proximal":
            ap.error("--out-dir is required for distal runs")
        args.out_dir = f"{M6A_RESULTS_ROOT}/m6A_switch_competence_{args.wl}_{args.wr}_guu_motif"
    return args


def main(argv=None):
    args = parse_args(argv)
    id_col, partner_col = id_columns(args.isoform)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"isoform={args.isoform} window=[{args.wl}:{args.wr}) source={args.source_dir} out_dir={args.out_dir}")

    print(f"Loading per_site.csv from the existing [{args.wl}:{args.wr}] run (reusing p(i,j) -- no RNAfold re-run) ...")
    per_site = pd.read_csv(os.path.join(args.source_dir, "per_site.csv"))
    print(f"  {len(per_site)} (transcript, m6A site) rows")

    needed_ids = set(per_site[id_col])
    pas_lookup = load_pas_lookup(GTF_PATH, needed_ids=needed_ids)
    anchor = anchor_of_run(args.source_dir)   # same window as the run whose folds are reused
    wl_of = window_wl_by_transcript(args.source_dir, id_col)   # per-transcript window starts, if any
    fasta = pysam.FastaFile(GENOME_FASTA)

    seq_cache = {}

    def get_seq(tx_id):
        if tx_id not in seq_cache:
            info = pas_lookup.loc[tx_id]
            seq, _ = fetch_window(fasta, info["chrom"], info["strand"], anchor_coordinate(info, anchor), wl_of.get(tx_id, args.wl), args.wr)
            seq_cache[tx_id] = ([""] + list(seq)) if seq is not None else None
        return seq_cache[tx_id]

    motif_type_col = []
    switch_competent_col = []
    n_done = 0

    for r in per_site.itertuples(index=False):
        n_done += 1
        if n_done % 5000 == 0:
            print(f"  processed {n_done} / {len(per_site)} sites ...")

        seq_1idx = get_seq(getattr(r, id_col))
        j = int(r.j)
        if seq_1idx is None or j <= 0:
            motif_type_col.append(np.nan)
            switch_competent_col.append(False)
            continue

        n = len(seq_1idx) - 1
        matches = motif_at_position(seq_1idx, j, n)
        switch_competent = bool(r.p_ij >= P_THRESHOLD and len(matches) > 0)
        motif_type_col.append(",".join(sorted(matches)) if matches else np.nan)
        switch_competent_col.append(switch_competent)

    per_site_out = per_site[[id_col, "gene", "chrom", "strand", "site_id", "i", "j", "p_ij",
                              "nuclear_lengthened"]].copy()
    for m in MOTIFS:
        per_site_out[f"motif_{m}"] = [
            (isinstance(mt, str) and m in mt.split(",")) for mt in motif_type_col
        ]
    per_site_out["motif_type"] = motif_type_col
    per_site_out["switch_competent"] = switch_competent_col

    print(f"\nswitch_competent breakdown: {per_site_out['switch_competent'].sum()} total / {len(per_site_out)}")
    print(per_site_out.loc[per_site_out["switch_competent"], "motif_type"].value_counts())

    per_tx_rows = []
    for tx_id, g in per_site_out.groupby(id_col):
        row0 = g.iloc[0]
        n_switch = {m: int((g[f"motif_{m}"] & g["switch_competent"]).sum()) for m in MOTIFS}
        n_switch_total = int(g["switch_competent"].sum())
        per_tx_rows.append({
            id_col: tx_id, "gene": row0["gene"], "chrom": row0["chrom"], "strand": row0["strand"],
            "n_m6a_sites": len(g),
            "n_switch_competent": n_switch_total,
            **{f"n_switch_competent_{m}": n_switch[m] for m in MOTIFS},
            "frac_switch_competent": (n_switch_total / len(g)) if len(g) > 0 else np.nan,
            "nuclear_lengthened": bool(row0["nuclear_lengthened"]),
        })
    per_tx_df = pd.DataFrame(per_tx_rows)

    # transcripts with zero m6A sites in window never appear in per_site.csv at all --
    # re-add them (n_m6a_sites=0) so per_transcript.csv covers all isoforms, matching the
    # convention of the main pipeline's per_transcript.csv (n_A carried over when the source
    # run has it -- legacy source runs predate that column)
    source_cols = pd.read_csv(os.path.join(args.source_dir, "per_transcript.csv"), nrows=0).columns
    usecols = [id_col, "gene", partner_col, "chrom", "strand", "nuclear_lengthened"]
    if "n_A" in source_cols:
        usecols.append("n_A")
    source_per_tx = pd.read_csv(os.path.join(args.source_dir, "per_transcript.csv"), usecols=usecols)[usecols]
    per_tx_df = source_per_tx.merge(per_tx_df.drop(columns=["gene", "chrom", "strand", "nuclear_lengthened"]),
                                     on=id_col, how="left")
    per_tx_df["n_m6a_sites"] = per_tx_df["n_m6a_sites"].fillna(0).astype(int)
    per_tx_df["n_switch_competent"] = per_tx_df["n_switch_competent"].fillna(0).astype(int)
    for m in MOTIFS:
        per_tx_df[f"n_switch_competent_{m}"] = per_tx_df[f"n_switch_competent_{m}"].fillna(0).astype(int)
    per_tx_df["frac_switch_competent"] = np.where(
        per_tx_df["n_m6a_sites"] > 0, per_tx_df["n_switch_competent"] / per_tx_df["n_m6a_sites"], np.nan,
    )

    per_site_path = os.path.join(args.out_dir, "per_site.csv")
    per_tx_path = os.path.join(args.out_dir, "per_transcript.csv")
    per_site_out.to_csv(per_site_path, index=False)
    per_tx_df.to_csv(per_tx_path, index=False)
    print(f"\nWrote {per_site_path} ({len(per_site_out)} rows)")
    print(f"Wrote {per_tx_path} ({len(per_tx_df)} rows)")

    write_params(
        args.out_dir,
        run_params={"step": "switch_competence_gu_motif", "isoform": args.isoform, "anchor": anchor, "wl": args.wl,
                    "wr": args.wr, "p_threshold": P_THRESHOLD, "motifs": MOTIFS,
                    "source_dir": os.path.abspath(args.source_dir)},
        script_paths=[os.path.abspath(__file__),
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), "competence.py")],
        input_paths={"gtf": GTF_PATH, "genome_fasta": GENOME_FASTA},
    )


if __name__ == "__main__":
    main()
