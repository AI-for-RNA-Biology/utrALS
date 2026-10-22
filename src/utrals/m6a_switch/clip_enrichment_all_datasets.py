#!/usr/bin/env python
"""
Run the switch-linked-element CLIP enrichment test for EVERY CLIP dataset in CLIP_ROOT
(eCLIP K562 + HepG2, iCLIP), in every region x element type, so the candidate RBPs can be placed
within the distribution of all RBPs.

Same test as clip_u_tract_switch_competent_enrichment.py / clip_guu_motif_switch_competent_enrichment.py,
reusing their outputs and functions so the numbers are identical for the candidates:
  - element catalogs: <root>/<region>/clip_u_tract_enrichment/u_tract_catalog.csv and
    <root>/<region>/clip_gu_motif_enrichment/motif_catalog.csv (written by the region pipeline)
  - null p0: switch-linked share of U-tract nucleotides (U-tracts) / of motif occurrences (GU motifs)
  - load_crosslinks / test_rbp from the U-tract script; one-sided binomial test
Each CLIP file is loaded once and tested in all region x element combinations. Region 0 (3'UTR
start) has no candidate list: its catalogs come from the CLIP scripts run with --catalog-only, and it
has no candidate rows / candidates_vs_rest test.

Datasets: every *.gff (not *.unmap) under CLIP_ROOT/eCLIP/{K652,HepG2} and CLIP_ROOT/iCLIP; when an
eCLIP dataset has both X.gff and X.sorted.gff, the unsorted one is used (as resolve_clip_file does).

Usage (a SLURM array splits the datasets into chunks, then one merge step):
  python clip_enrichment_all_datasets.py --list                      # write datasets.csv
  python clip_enrichment_all_datasets.py --chunk K --n-chunks N      # test chunk K -> parts/part_K.csv
  python clip_enrichment_all_datasets.py --merge                     # -> clip_enrichment_all_datasets.csv,
                                                                     #    candidates_vs_rest.csv
fold_change = observed_rate / expected_rate (0 when no crosslink falls in a switch-linked element).
candidates_vs_rest.csv: per region x element, do the region's candidate datasets have higher
fold changes than all other datasets? (one-sided Mann-Whitney; median FC and median percentile
of the candidates). Same for -log10(binomial p) (columns *_neglog10p; p-values below 1e-300 floored).
Output dir: <root>/clip_all_datasets/
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .competence import M6A_RESULTS_ROOT, write_params
from .clip_enrichment_u_tract import CLIP_ROOT, load_crosslinks, test_rbp

DEFAULT_ROOT = os.path.join(M6A_RESULTS_ROOT, "regions")
REGIONS = {"3'UTR start [-100:+100]": "utr_start_m100_100", "proximal [-350:0]": "proximal_m350_0", "proximal [0:200]": "proximal_0_200",
           "distal [-350:0]": "distal_m350_0"}
# element -> (clip subdir, catalog file, start col, end col, p0 weighting)
ELEMENTS = {
    "U-tracts": ("clip_u_tract_enrichment", "u_tract_catalog.csv", "u_run_genomic_start", "u_run_genomic_end", "nt"),
    "GU motifs": ("clip_gu_motif_enrichment", "motif_catalog.csv", "occurrence_genomic_start", "occurrence_genomic_end", "count"),
}


def list_datasets():
    rows = []
    for source in ("K652", "HepG2"):
        files = sorted(glob.glob(os.path.join(CLIP_ROOT, "eCLIP", source, "*.gff")))
        unsorted = {f for f in files if not f.endswith(".sorted.gff")}
        for f in files:
            if f.endswith(".sorted.gff") and f.replace(".sorted.gff", ".gff") in unsorted:
                continue
            name = os.path.basename(f)
            prefix = name.split("_ENC")[0]                       # e.g. "DDX6_2"
            rows.append({"dataset": f"{prefix}_{source}", "RBP": prefix.rsplit("_", 1)[0].upper(),
                         "assay": f"eCLIP {source}", "path": f})
    for f in sorted(glob.glob(os.path.join(CLIP_ROOT, "iCLIP", "*.gff"))):
        name = os.path.basename(f)
        label = re.sub(r"^group_\d+_", "", name).split("_sum")[0].removesuffix("hg38").rstrip("-")
        rbp = re.sub(r"^all-", "", label).split("-")[0].upper()
        rows.append({"dataset": f"{label}_iCLIP", "RBP": rbp, "assay": "iCLIP", "path": f})
    df = pd.DataFrame(rows)
    df["size_bytes"] = [os.path.getsize(p) for p in df["path"]]
    assert df["dataset"].is_unique and df["path"].is_unique
    return df


def load_targets(root):
    """[(region, element, p0, {(chrom, strand): (starts, ends, categories)})]"""
    targets = []
    for region, sub in REGIONS.items():
        for element, (cdir, cat_file, s_col, e_col, weight) in ELEMENTS.items():
            cat = pd.read_csv(os.path.join(root, sub, cdir, cat_file))
            w = (cat[e_col] - cat[s_col]) if weight == "nt" else pd.Series(1, index=cat.index)
            p0 = w[cat["category"] == "switch"].sum() / w.sum()
            grouped = {}
            for (chrom, strand), g in cat.groupby(["chrom", "strand"]):
                g = g.sort_values(s_col)
                grouped[(chrom, strand)] = (g[s_col].values, g[e_col].values, g["category"].values)
            targets.append((region, element, p0, grouped))
    return targets


def run_chunk(root, out_dir, k, n_chunks):
    ds = pd.read_csv(os.path.join(out_dir, "datasets.csv"))
    # balance chunks by file size: largest first, round-robin
    order = ds.sort_values("size_bytes", ascending=False).index
    mine = ds.loc[[i for j, i in enumerate(order) if j % n_chunks == k]]
    targets = load_targets(root)
    rows = []
    for d in mine.itertuples(index=False):
        print(f"[{k}] {d.dataset} ({d.size_bytes / 1e6:.0f} MB)", flush=True)
        clip_by_key = load_crosslinks(d.path)
        for region, element, p0, grouped in targets:
            n, kk = test_rbp(clip_by_key, grouped)
            obs = kk / n if n else np.nan
            rows.append({"dataset": d.dataset, "RBP": d.RBP, "assay": d.assay, "region": region,
                         "element": element, "n_crosslinks_in_elements": n, "n_in_switch_elements": kk,
                         "observed_rate": obs, "expected_rate": p0,
                         "odds_ratio": (obs / (1 - obs)) / (p0 / (1 - p0)) if 0 < obs < 1 else np.nan,
                         "p_value": binomtest(kk, n, p0, alternative="greater").pvalue if n else np.nan})
        del clip_by_key
    os.makedirs(os.path.join(out_dir, "parts"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "parts", f"part_{k:03d}.csv"), index=False)


def merge(root, out_dir, n_chunks):
    parts = sorted(glob.glob(os.path.join(out_dir, "parts", "part_*.csv")))
    if len(parts) != n_chunks:
        sys.exit(f"ERROR: {len(parts)}/{n_chunks} parts present")
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    ds = pd.read_csv(os.path.join(out_dir, "datasets.csv"))
    assert set(df["dataset"]) == set(ds["dataset"]) and len(df) == len(ds) * len(REGIONS) * len(ELEMENTS)
    df["fold_change"] = df["observed_rate"] / df["expected_rate"]
    df["log2_odds_ratio"] = np.log2(df["odds_ratio"])
    # candidate flag: the region's own candidate list (read from its CLIP run), matched by file path
    from .clip_enrichment_u_tract import resolve_clip_file
    path_of = dict(zip(ds["dataset"], ds["path"]))
    df["is_candidate"] = False
    for region, sub in REGIONS.items():
        cand_file = os.path.join(root, sub, "clip_u_tract_enrichment", "clip_enrichment.csv")
        if not os.path.exists(cand_file):   # catalog-only region: no candidates
            continue
        cand = pd.read_csv(cand_file)
        cand_paths = {resolve_clip_file(c) for c in cand["CLIP"]}
        m = df["region"] == region
        df.loc[m, "is_candidate"] = df.loc[m, "dataset"].map(path_of).isin(cand_paths)
    for (region, element), g in df.groupby(["region", "element"]):
        ok = g["p_value"].notna()
        df.loc[g.index[ok], "q_value_bh"] = multipletests(g.loc[ok, "p_value"], method="fdr_bh")[1]
        df.loc[g.index, "rank_fc"] = g["fold_change"].rank(ascending=False, method="min")
        df.loc[g.index, "percentile_fc"] = 100 * g["fold_change"].rank(pct=True)
    df = df.sort_values(["region", "element", "rank_fc"])
    out = os.path.join(out_dir, "clip_enrichment_all_datasets.csv")
    df.to_csv(out, index=False)
    print(f"Wrote {out} ({len(df)} rows, {df['dataset'].nunique()} datasets)")

    rows = []
    for (region, element), g in df.groupby(["region", "element"], sort=False):
        g = g[g["fold_change"].notna()]
        c, r = g[g["is_candidate"]], g[~g["is_candidate"]]
        if c.empty:   # region without candidate list
            continue
        u, p = mannwhitneyu(c["fold_change"], r["fold_change"], alternative="greater")
        nlp_c, nlp_r = (-np.log10(x["p_value"].clip(lower=1e-300)) for x in (c, r))
        u2, p2 = mannwhitneyu(nlp_c, nlp_r, alternative="greater")
        rows.append({"region": region, "element": element, "n_candidates": len(c), "n_rest": len(r),
                     "median_fc_candidates": c["fold_change"].median(),
                     "median_fc_rest": r["fold_change"].median(),
                     "median_percentile_candidates": c["percentile_fc"].median(),
                     "frac_rest_significant_q05": (r["q_value_bh"] < 0.05).mean(),
                     "frac_candidates_significant_q05": (c["q_value_bh"] < 0.05).mean(),
                     "test": "Mann-Whitney U (candidates > rest)", "statistic": u, "p_value": p,
                     "median_neglog10p_candidates": nlp_c.median(), "median_neglog10p_rest": nlp_r.median(),
                     "statistic_neglog10p": u2, "p_value_neglog10p": p2})
    summ = pd.DataFrame(rows)
    for pcol, qcol in [("p_value", "q_value_bh"), ("p_value_neglog10p", "q_value_bh_neglog10p")]:
        summ[qcol] = multipletests(summ[pcol], method="fdr_bh")[1]   # BH over the region x element tests
    summ.to_csv(os.path.join(out_dir, "candidates_vs_rest.csv"), index=False)
    print(summ.to_string(index=False))
    write_params(out_dir, run_params={"step": "clip_enrichment_all_datasets", "n_chunks": n_chunks,
                                      "regions": REGIONS, "elements": list(ELEMENTS)},
                 script_paths=[os.path.abspath(__file__)] + [
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), s)
                     for s in ("clip_enrichment_u_tract.py", "competence.py")],
                 input_paths={"root": os.path.abspath(root), "clip_root": CLIP_ROOT})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--n-chunks", type=int, default=40)
    ap.add_argument("--merge", action="store_true")
    args = ap.parse_args(argv)
    out_dir = os.path.join(args.root, "clip_all_datasets")
    os.makedirs(out_dir, exist_ok=True)
    if args.list:
        ds = list_datasets()
        ds.to_csv(os.path.join(out_dir, "datasets.csv"), index=False)
        print(f"{len(ds)} datasets: " + ds["assay"].value_counts().to_string().replace("\n", ", "))
    elif args.chunk is not None:
        run_chunk(args.root, out_dir, args.chunk, args.n_chunks)
    elif args.merge:
        merge(args.root, out_dir, args.n_chunks)
    else:
        ap.error("one of --list, --chunk, --merge")


if __name__ == "__main__":
    main()
