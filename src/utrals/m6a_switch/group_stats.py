#!/usr/bin/env python
"""
Nuclear-lengthened vs. other statistics for one switch-competence run directory (U-tract run
from m6a_switch_competence_analysis.py, or GU-motif run from
m6a_switch_competence_guu_motif_analysis.py), so notebooks only have to read and plot.

Requires per_transcript.csv to have n_A (runs produced by the parameterized pipeline do;
legacy directories predating it don't).

Writes into --run-dir:
  - group_stats.csv            : one row per (motif, metric) comparison. motif is "all" for
                                  the run's own switch-competence call, or UGU/UUG/GUG/GUU for
                                  the per-motif breakdown of a GU-motif run (not mutually
                                  exclusive). Metrics:
      raw_per_transcript         mean # switch-competent sites per transcript, all
                                  transcripts; two-sided Mann-Whitney
      frac_per_m6a_site          mean per-transcript fraction of m6A sites that are
                                  switch-competent, transcripts with >=1 m6A site; Mann-Whitney
      per_A_per_transcript       mean per-transcript switch-competent sites per A; Mann-Whitney
      pooled_per_m6a_site        pooled switch-competent / m6A sites; Fisher's exact
      pooled_per_A               pooled switch-competent / A nucleotides; Fisher's exact
    Odds ratios are nuclear-lengthened vs. other (OR > 1 => enriched in nuclear-lengthened).
    q_value_bh: Benjamini-Hochberg over the 5 metrics of each motif (alternative normalizations
    of the same nuclear-lengthened vs. other comparison).
  - per_transcript_summary.csv : per (motif, group) n, mean, SEM and 95% CI half-width
                                  (1.96 x SEM) of switch-competent sites per transcript, for
                                  bar plots with error bars.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .competence import read_run_params, write_params

GROUPS = {True: "Nuclear-lengthened", False: "Other"}


def _mwu(values, is_nl):
    values, is_nl = np.asarray(values, dtype=float), np.asarray(is_nl, dtype=bool)
    keep = ~np.isnan(values)
    values, is_nl = values[keep], is_nl[keep]
    u, p = mannwhitneyu(values[is_nl], values[~is_nl], alternative="two-sided")
    return dict(n_nuclear_lengthened=int(is_nl.sum()), n_other=int((~is_nl).sum()),
                value_nuclear_lengthened=values[is_nl].mean(), value_other=values[~is_nl].mean(),
                odds_ratio=np.nan, statistic=u, p_value=p, test="Mann-Whitney U (two-sided)")


def _fisher(k_nl, n_nl, k_other, n_other):
    """k successes out of n trials per group; OR = odds(NL) / odds(Other)."""
    table = [[k_nl, n_nl - k_nl], [k_other, n_other - k_other]]
    odds, p = fisher_exact(table, alternative="two-sided")
    return dict(n_nuclear_lengthened=int(n_nl), n_other=int(n_other),
                value_nuclear_lengthened=k_nl / n_nl, value_other=k_other / n_other,
                odds_ratio=odds, statistic=np.nan, p_value=p, test="Fisher's exact (two-sided)")


def group_stats(per_site, per_tx, count_col, site_mask, motif):
    nl_tx = per_tx["nuclear_lengthened"].astype(bool)
    nl_site = per_site["nuclear_lengthened"].astype(bool)
    has_sites = per_tx["n_m6a_sites"] > 0
    counts = per_tx[count_col]

    rows = [
        dict(metric="raw_per_transcript", **_mwu(counts, nl_tx)),
        dict(metric="frac_per_m6a_site", **_mwu((counts / per_tx["n_m6a_sites"])[has_sites], nl_tx[has_sites])),
        dict(metric="per_A_per_transcript", **_mwu(counts / per_tx["n_A"].replace(0, np.nan), nl_tx)),
        dict(metric="pooled_per_m6a_site", **_fisher(
            int(site_mask[nl_site].sum()), int(nl_site.sum()),
            int(site_mask[~nl_site].sum()), int((~nl_site).sum()))),
        dict(metric="pooled_per_A", **_fisher(
            int(counts[nl_tx].sum()), int(per_tx.loc[nl_tx, "n_A"].sum()),
            int(counts[~nl_tx].sum()), int(per_tx.loc[~nl_tx, "n_A"].sum()))),
    ]
    for r in rows:
        r["motif"] = motif
    return rows


def per_transcript_summary(per_tx, count_col, motif):
    rows = []
    for is_nl, label in GROUPS.items():
        v = per_tx.loc[per_tx["nuclear_lengthened"].astype(bool) == is_nl, count_col]
        sem = v.sem()
        rows.append(dict(motif=motif, group=label, n=len(v), mean=v.mean(), sem=sem, ci95=1.96 * sem))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="switch-competence run directory")
    args = ap.parse_args(argv)

    per_site = pd.read_csv(os.path.join(args.run_dir, "per_site.csv"))
    per_tx = pd.read_csv(os.path.join(args.run_dir, "per_transcript.csv"))
    if "n_A" not in per_tx.columns:
        sys.exit(f"{args.run_dir}/per_transcript.csv has no n_A column -- re-run it with the "
                 "parameterized pipeline (legacy runs predate n_A)")

    motifs = [c[len("motif_"):] for c in per_site.columns if c.startswith("motif_") and c != "motif_type"]
    comparisons = [("all", "n_switch_competent", per_site["switch_competent"].astype(bool))]
    for m in motifs:
        comparisons.append((m, f"n_switch_competent_{m}",
                            per_site["switch_competent"].astype(bool) & per_site[f"motif_{m}"].astype(bool)))

    stats_rows, summary_rows = [], []
    for motif, count_col, site_mask in comparisons:
        stats_rows += group_stats(per_site, per_tx, count_col, site_mask, motif)
        summary_rows += per_transcript_summary(per_tx, count_col, motif)

    cols = ["motif", "metric", "test", "n_nuclear_lengthened", "n_other", "value_nuclear_lengthened",
            "value_other", "odds_ratio", "statistic", "p_value"]
    stats_df = pd.DataFrame(stats_rows)[cols]
    stats_df["q_value_bh"] = np.nan
    for _, g in stats_df.groupby("motif", sort=False):
        ok = g["p_value"].notna()
        stats_df.loc[g.index[ok], "q_value_bh"] = multipletests(g.loc[ok, "p_value"], method="fdr_bh")[1]
    summary_df = pd.DataFrame(summary_rows)
    stats_df.to_csv(os.path.join(args.run_dir, "group_stats.csv"), index=False)
    summary_df.to_csv(os.path.join(args.run_dir, "per_transcript_summary.csv"), index=False)
    print(stats_df.to_string(index=False))
    print(summary_df.to_string(index=False))

    # record that stats were computed, without overwriting the run's own params.json
    run = read_run_params(args.run_dir) or {}
    stats_dir = os.path.join(args.run_dir, "stats_provenance")
    os.makedirs(stats_dir, exist_ok=True)
    write_params(stats_dir, run_params={"step": "group_stats", "run_dir": os.path.abspath(args.run_dir),
                                        "source_run": run},
                 script_paths=[os.path.abspath(__file__)], input_paths={})


if __name__ == "__main__":
    main()
