#!/usr/bin/env python
"""
Collect and compare the m6A switch-competence pipeline outputs of several regions (one
run_m6a_switch_pipeline.sh output root per region) into <root>/comparison/:

  overview.csv                 per region x definition (UUUs / UG motifs) x group (all,
                               Nuclear-lengthened, Other): transcripts, transcripts with >= 1
                               m6A site, m6A sites, A nucleotides, switch-competent sites, and
                               the pooled rates per m6A site and per 1,000 A
  nl_vs_other.csv              every region's group_stats.csv (nuclear-lengthened vs. other),
                               stacked with region/definition columns
  nl_vs_other_pooled_or.csv    pooled odds ratio (nuclear-lengthened vs. other) of a site being
                               switch-competent, per m6A site and per A nucleotide, with Woolf 95% CI
                               and Fisher p (same tests as group_stats.csv, plus the CI)
  per_transcript_summary.csv   every region's per_transcript_summary.csv, stacked
  per_transcript_frac_summary.csv  mean per-transcript fraction of m6A sites that are
                               switch-competent (transcripts with >= 1 m6A site; the values the
                               frac_per_m6a_site Mann-Whitney test compares), with SEM and 95% CI
  vcp_d35_per_transcript_summary.csv  pairs split by the D35 nuclear VCP vs. CTRL shift
                               (VCP-shortened = pairs_proximals_D35_Nuclear_Wildtype_vs_D35_Nuclear_VCP.csv,
                               VCP-lengthened = pairs_distals_..., Other = rest of the analysed pairs):
                               mean +/- SEM / 95% CI of switch-competent sites per transcript
                               (metric raw_per_transcript, all transcripts) and of the fraction of
                               m6A sites that are switch-competent (frac_per_m6a_site, transcripts
                               with >= 1 m6A site)
  vcp_d35_pairwise_tests.csv   the three pairwise two-sided Mann-Whitney tests between those groups,
                               per region x definition x metric, with BH q over the 3 tests
  nuc_d22_d35_per_transcript_summary.csv, nuc_d22_d35_pairwise_tests.csv
                               the same for the D22 vs. D35 nuclear wild-type shift:
                               Nuclear-shortened (pairs_proximals_D22_Nuclear_Wildtype_vs_D35_Nuclear_Wildtype.csv),
                               Nuclear-lengthened (pairs_distals_..., the nuclear_lengthened label of
                               every other table), Other = rest of the analysed pairs
  nl_vcp_shortened_per_transcript_summary.csv, nl_vcp_shortened_tests.csv
                               among the nuclear-lengthened pairs (CTRL D22 vs. D35 nuclear) only:
                               those also shortened in VCP vs. CTRL (D35 nuclear) vs. the rest;
                               same two metrics, two-sided Mann-Whitney, BH over all region x
                               definition x metric tests
  between_region_tests.csv     all transcripts pooled: is the switch-competent rate different
                               between regions? Chi-square across all regions and pairwise
                               Fisher tests (BH-corrected within definition x metric), per m6A
                               site and per A nucleotide
  clip_enrichment_all_regions.csv  every region's CLIP enrichment (U-tracts and GU motifs),
                               stacked, with log2 OR and BH q-values within region x element type

Regions are defined in REGIONS below (label -> subdirectory of --root); each region's window,
isoform and candidate list are read back from its params.json files.

Usage: python compare_switch_regions.py [--root <dir>]
"""

import argparse
import itertools
import json
import os
import sys
from ..paths import PROJECT_DIR as _PROJECT_DIR  # data locations, see paths.env

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .competence import M6A_RESULTS_ROOT, write_params

DEFAULT_ROOT = os.path.join(M6A_RESULTS_ROOT, "regions")
# label -> subdirectory; order = upstream-to-downstream along the 3'UTR
REGIONS = {
    "3'UTR start [-100:+100]": "utr_start_m100_100",   # R0: around the stop codon (no CLIP candidates)
    "proximal [-350:0]": "proximal_m350_0",
    "proximal [0:200]": "proximal_0_200",
    "distal [-350:0]": "distal_m350_0",
}
SHIFTS_DIR = f"{_PROJECT_DIR}/MN_diff/results/3UTR_analysis/shifts"
# three-way pair groupings: name -> {group: shift file}; pairs in neither file are "Other"
GROUPINGS = {
    "vcp_d35": {
        "VCP-shortened": os.path.join(SHIFTS_DIR, "pairs_proximals_D35_Nuclear_Wildtype_vs_D35_Nuclear_VCP.csv"),
        "VCP-lengthened": os.path.join(SHIFTS_DIR, "pairs_distals_D35_Nuclear_Wildtype_vs_D35_Nuclear_VCP.csv"),
    },
    "nuc_d22_d35": {
        "Nuclear-shortened": os.path.join(SHIFTS_DIR, "pairs_proximals_D22_Nuclear_Wildtype_vs_D35_Nuclear_Wildtype.csv"),
        "Nuclear-lengthened": os.path.join(SHIFTS_DIR, "pairs_distals_D22_Nuclear_Wildtype_vs_D35_Nuclear_Wildtype.csv"),
    },
}
DEFINITIONS = {"UUUs": "switch_competence_u_tract", "UG motifs": "switch_competence_gu_motif"}
CLIP = {"U-tracts": ("clip_u_tract_enrichment", "n_crosslinks_in_u_tracts", "n_in_switch_u_tract"),
        "GU motifs": ("clip_gu_motif_enrichment", "n_crosslinks_in_motifs", "n_in_switch_motif")}


def _params(run_dir):
    with open(os.path.join(run_dir, "params.json")) as f:
        return json.load(f)


def overview_rows(region, definition, per_tx):
    rows = []
    groups = {"all": per_tx,
              "Nuclear-lengthened": per_tx[per_tx["nuclear_lengthened"]],
              "Other": per_tx[~per_tx["nuclear_lengthened"]]}
    for g, d in groups.items():
        n_sites, n_sw, n_a = int(d["n_m6a_sites"].sum()), int(d["n_switch_competent"].sum()), int(d["n_A"].sum())
        rows.append({
            "region": region, "definition": definition, "group": g,
            "n_transcripts": len(d), "n_transcripts_with_m6a": int((d["n_m6a_sites"] > 0).sum()),
            "n_m6a_sites": n_sites, "n_A": n_a, "n_switch_competent": n_sw,
            "frac_switch_per_m6a_site": n_sw / n_sites if n_sites else np.nan,
            "switch_per_1000_A": 1000 * n_sw / n_a if n_a else np.nan,
            "mean_switch_per_transcript": d["n_switch_competent"].mean(),
        })
    return rows


def frac_summary_rows(region, definition, per_tx):
    d = per_tx[per_tx["n_m6a_sites"] > 0]
    frac = d["n_switch_competent"] / d["n_m6a_sites"]
    rows = []
    for g, m in [("Nuclear-lengthened", d["nuclear_lengthened"]), ("Other", ~d["nuclear_lengthened"])]:
        v = frac[m]
        rows.append({"region": region, "definition": definition, "group": g, "n": len(v),
                     "mean": v.mean(), "sem": v.sem(), "ci95": 1.96 * v.sem()})
    return rows


def load_group_pairs(shift_files):
    """{(proximal_id, distal_id): group} for the pairs listed in each group's shift file."""
    labels = {}
    for group, path in shift_files.items():
        for pair in pd.read_csv(path, usecols=["proximal_id", "distal_id"]).itertuples(index=False):
            if tuple(pair) in labels:
                raise ValueError(f"pair {tuple(pair)} is in more than one shift file of {list(shift_files)}")
            labels[tuple(pair)] = group
    return labels


def three_group_rows(region, definition, per_tx, pair_groups, groups):
    """Per-group summaries and pairwise Mann-Whitney tests for a three-way pair grouping."""
    group = pd.Series([pair_groups.get(k, "Other") for k in zip(per_tx["proximal_id"], per_tx["distal_id"])],
                      index=per_tx.index)
    has_sites = per_tx["n_m6a_sites"] > 0
    metrics = {
        "raw_per_transcript": per_tx["n_switch_competent"],
        "frac_per_m6a_site": (per_tx["n_switch_competent"] / per_tx["n_m6a_sites"])[has_sites],
    }
    summ, tests = [], []
    for metric, values in metrics.items():
        g = group.reindex(values.index)
        for name in groups:
            v = values[g == name]
            summ.append({"region": region, "definition": definition, "metric": metric, "group": name,
                         "n": len(v), "mean": v.mean(), "sem": v.sem(), "ci95": 1.96 * v.sem()})
        rows = []
        for a, b in itertools.combinations(groups, 2):
            va, vb = values[g == a], values[g == b]
            u, p = mannwhitneyu(va, vb, alternative="two-sided")
            rows.append({"region": region, "definition": definition, "metric": metric,
                         "group_a": a, "group_b": b, "n_a": len(va), "n_b": len(vb),
                         "mean_a": va.mean(), "mean_b": vb.mean(),
                         "test": "Mann-Whitney U (two-sided)", "statistic": u, "p_value": p})
        q = multipletests([r["p_value"] for r in rows], method="fdr_bh")[1]
        for r, qi in zip(rows, q):
            r["q_value_bh"] = qi
        tests += rows
    return summ, tests


def nl_vcp_shortened_rows(region, definition, per_tx, vcp_short_pairs):
    """Among nuclear-lengthened pairs: VCP-shortened (D35 nuclear) vs. not."""
    d = per_tx[per_tx["nuclear_lengthened"]]
    is_short = pd.Series([k in vcp_short_pairs for k in zip(d["proximal_id"], d["distal_id"])], index=d.index)
    has_sites = d["n_m6a_sites"] > 0
    metrics = {"raw_per_transcript": d["n_switch_competent"],
               "frac_per_m6a_site": (d["n_switch_competent"] / d["n_m6a_sites"])[has_sites]}
    summ, tests = [], []
    for metric, v in metrics.items():
        s = is_short.reindex(v.index)
        groups = {"NL, VCP-shortened": v[s], "NL, not VCP-shortened": v[~s]}
        for g, x in groups.items():
            summ.append({"region": region, "definition": definition, "metric": metric, "group": g,
                         "n": len(x), "mean": x.mean(), "sem": x.sem(), "ci95": 1.96 * x.sem()})
        a, b = groups.values()
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        tests.append({"region": region, "definition": definition, "metric": metric,
                      "group_a": "NL, VCP-shortened", "group_b": "NL, not VCP-shortened", "n_a": len(a), "n_b": len(b),
                      "mean_a": a.mean(), "mean_b": b.mean(), "test": "Mann-Whitney U (two-sided)",
                      "statistic": u, "p_value": p})
    return summ, tests


def nl_vs_other_or(overview):
    rows = []
    for (region, definition), d in overview.groupby(["region", "definition"], sort=False):
        d = d.set_index("group")
        for metric, denom in [("per_m6a_site", "n_m6a_sites"), ("per_A", "n_A")]:
            a = d.loc["Nuclear-lengthened", "n_switch_competent"]; b = d.loc["Nuclear-lengthened", denom] - a
            c = d.loc["Other", "n_switch_competent"]; e = d.loc["Other", denom] - c
            odds, p = fisher_exact([[a, b], [c, e]])
            se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / e) if min(a, b, c, e) > 0 else np.nan
            rows.append({"region": region, "definition": definition, "metric": metric,
                         "rate_nuclear_lengthened": a / (a + b), "rate_other": c / (c + e),
                         "odds_ratio": odds, "ci95_low": np.exp(np.log(odds) - 1.96 * se),
                         "ci95_high": np.exp(np.log(odds) + 1.96 * se), "p_value": p})
    out = pd.DataFrame(rows)
    # q: the same Fisher tests are in each run's group_stats.csv (pooled_per_m6a_site / pooled_per_A);
    # filled in main() from there so both tables carry the same q
    return out


def add_pooled_q(orr, nl):
    """q of each pooled OR = q of the identical Fisher test in group_stats.csv (motif 'all')."""
    q = nl[nl["motif"] == "all"].set_index(["region", "definition", "metric"])["q_value_bh"]
    keys = zip(orr["region"], orr["definition"], "pooled_" + orr["metric"])
    orr["q_value_bh"] = [q.get(k, np.nan) for k in keys]
    return orr


def between_region_tests(overview):
    rows = []
    allg = overview[overview["group"] == "all"]
    for definition, d in allg.groupby("definition", sort=False):
        d = d.set_index("region")
        for metric, denom in [("per_m6a_site", "n_m6a_sites"), ("per_A", "n_A")]:
            table = np.array([[d.loc[r, "n_switch_competent"], d.loc[r, denom] - d.loc[r, "n_switch_competent"]]
                              for r in d.index])
            chi2, p, dof, _ = chi2_contingency(table)
            rows.append({"definition": definition, "metric": metric, "comparison": "all regions",
                         "test": f"chi-square (dof={dof})", "rate_a": np.nan, "rate_b": np.nan,
                         "odds_ratio": np.nan, "statistic": chi2, "p_value": p})
            pair_rows = []
            for a, b in itertools.combinations(d.index, 2):
                ka, na = d.loc[a, "n_switch_competent"], d.loc[a, denom]
                kb, nb = d.loc[b, "n_switch_competent"], d.loc[b, denom]
                odds, p = fisher_exact([[ka, na - ka], [kb, nb - kb]])
                pair_rows.append({"definition": definition, "metric": metric, "comparison": f"{a} vs {b}",
                                  "test": "Fisher exact (two-sided)", "rate_a": ka / na, "rate_b": kb / nb,
                                  "odds_ratio": odds, "statistic": np.nan, "p_value": p})
            q = multipletests([r["p_value"] for r in pair_rows], method="fdr_bh")[1]
            for r, qi in zip(pair_rows, q):
                r["q_value_bh"] = qi
            rows += pair_rows
    out = pd.DataFrame(rows)
    chi = out["comparison"] == "all regions"   # BH over the chi-square tests (definition x metric)
    out.loc[chi, "q_value_bh"] = multipletests(out.loc[chi, "p_value"], method="fdr_bh")[1]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT, help=f"directory holding the region subdirectories (default {DEFAULT_ROOT})")
    args = ap.parse_args(argv)
    out_dir = os.path.join(args.root, "comparison")
    os.makedirs(out_dir, exist_ok=True)

    overview, nl, summ, frac_summ, clip, region_meta = [], [], [], [], [], []
    grouping_pairs = {k: load_group_pairs(files) for k, files in GROUPINGS.items()}
    grouping_summ = {k: [] for k in GROUPINGS}
    grouping_tests = {k: [] for k in GROUPINGS}
    vcp_short_pairs = set(load_group_pairs({"VCP-shortened": GROUPINGS["vcp_d35"]["VCP-shortened"]}))
    nlvs_summ, nlvs_tests = [], []
    for region, sub in REGIONS.items():
        rdir = os.path.join(args.root, sub)
        run = _params(os.path.join(rdir, DEFINITIONS["UUUs"]))["run"]
        has_clip = os.path.exists(os.path.join(rdir, CLIP["U-tracts"][0], "clip_enrichment.csv"))   # candidate test run
        cand = _params(os.path.join(rdir, CLIP["U-tracts"][0]))["inputs"]["candidates"] if has_clip else None
        region_meta.append({"region": region, "subdir": sub, "isoform": run["isoform"],
                            "anchor": run.get("anchor", "pas"), "wl": run["wl"], "wr": run["wr"],
                            "n_windows_cut_at_exon_start": run.get("n_windows_cut", 0), "candidates": cand})
        for definition, ddir in DEFINITIONS.items():
            run_dir = os.path.join(rdir, ddir)
            per_tx = pd.read_csv(os.path.join(run_dir, "per_transcript.csv"))
            per_tx["nuclear_lengthened"] = per_tx["nuclear_lengthened"].astype(bool)
            overview += overview_rows(region, definition, per_tx)
            frac_summ += frac_summary_rows(region, definition, per_tx)
            for k, files in GROUPINGS.items():
                s_rows, t_rows = three_group_rows(region, definition, per_tx, grouping_pairs[k],
                                                  ["Other"] + list(files))
                grouping_summ[k] += s_rows
                grouping_tests[k] += t_rows
            s_rows, t_rows = nl_vcp_shortened_rows(region, definition, per_tx, vcp_short_pairs)
            nlvs_summ += s_rows
            nlvs_tests += t_rows
            for name, acc in [("group_stats.csv", nl), ("per_transcript_summary.csv", summ)]:
                df = pd.read_csv(os.path.join(run_dir, name))
                df.insert(0, "definition", definition)
                df.insert(0, "region", region)
                acc.append(df)
        for element, (cdir, n_col, k_col) in (CLIP.items() if has_clip else []):
            df = pd.read_csv(os.path.join(rdir, cdir, "clip_enrichment.csv"))
            df = df.rename(columns={n_col: "n_crosslinks_in_elements", k_col: "n_in_switch_elements"})
            df.insert(0, "element", element)
            df.insert(0, "region", region)
            df["candidates"] = os.path.basename(cand)
            df["log2_odds_ratio"] = np.log2(df["odds_ratio"])
            ok = df["p_value"].notna()
            df.loc[ok, "q_value_bh"] = multipletests(df.loc[ok, "p_value"], method="fdr_bh")[1]
            clip.append(df)

    overview = pd.DataFrame(overview)
    outputs = {
        "regions.csv": pd.DataFrame(region_meta),
        "overview.csv": overview,
        "nl_vs_other.csv": pd.concat(nl, ignore_index=True),
        "nl_vs_other_pooled_or.csv": add_pooled_q(nl_vs_other_or(overview), pd.concat(nl, ignore_index=True)),
        "per_transcript_summary.csv": pd.concat(summ, ignore_index=True),
        "per_transcript_frac_summary.csv": pd.DataFrame(frac_summ),
        **{f"{k}_per_transcript_summary.csv": pd.DataFrame(v) for k, v in grouping_summ.items()},
        **{f"{k}_pairwise_tests.csv": pd.DataFrame(v) for k, v in grouping_tests.items()},
        "nl_vcp_shortened_per_transcript_summary.csv": pd.DataFrame(nlvs_summ),
        "nl_vcp_shortened_tests.csv": pd.DataFrame(nlvs_tests).assign(
            q_value_bh=lambda t: multipletests(t["p_value"], method="fdr_bh")[1]),
        "between_region_tests.csv": between_region_tests(overview),
        "clip_enrichment_all_regions.csv": pd.concat(clip, ignore_index=True),
    }
    for name, df in outputs.items():
        df.to_csv(os.path.join(out_dir, name), index=False)
        print(f"Wrote {os.path.join(out_dir, name)} ({len(df)} rows)")
    print(overview[overview["group"] == "all"].to_string(index=False))

    write_params(
        out_dir,
        run_params={"step": "compare_switch_regions", "regions": REGIONS},
        script_paths=[os.path.abspath(__file__)],
        input_paths={"root": os.path.abspath(args.root),
                     **{f"{k}:{g}": f for k, files in GROUPINGS.items() for g, f in files.items()}},
    )


if __name__ == "__main__":
    main()
