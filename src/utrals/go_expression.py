"""Expression of a GO biological-process gene set over the differentiation time course
(D0-D35), one compartment at a time, CTRL (WT) and VCP handled separately.

Gene set: human (taxon 9606) "involved_in" annotations of the GO term *and its is_a /
part_of descendants*, fetched from QuickGO (all evidence codes, incl. IEA). Unreviewed
TrEMBL entries without an HGNC symbol (symbol == UniProt accession / DKFZp clone) are
dropped, and outdated aliases are mapped to the current symbol.

Per-gene values come from utrals.expression.gene_expression_over_time(), i.e.
log2FC vs. that same replicate's own D0 value, plus plain log2(expression+1). Two pathway
scores per sample:
- median log2FC vs. D0 across the genes of the set;
- mean per-gene z-score of log2(expression+1), each gene z-scored across all samples of the
  chosen compartment (both genotypes), so CTRL and VCP live on the same scale.
Trend over time: Spearman correlation of the value with Day, over the 4 replicates x 6 days
of one genotype (n = 24), BH-corrected across genes.
"""

import io
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import spearmanr, t as t_dist, ttest_ind
from statsmodels.stats.multitest import multipletests

from .expression import DAY_ORDER, GENOTYPE_PALETTE, gene_expression_over_time

GENOTYPE_ORDER = ["WT", "VCP"]
GENOTYPE_LABEL = {"WT": "CTRL", "VCP": "VCP"}
SYMBOL_ALIASES = {"SHFM1": "SEM1"}

QUICKGO_URL = (
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
    "?goId={go_id}&goUsage=descendants&goUsageRelationships=is_a,part_of"
    "&taxonId=9606&taxonUsage=exact&geneProductType=protein"
    "&selectedFields=geneProductId,symbol,qualifier,goId,goName,evidenceCode,goEvidence,reference,assignedBy"
    "&downloadLimit=50000"
)


def fetch_go_annotations(go_id, cache_path=None, timeout=300):
    """QuickGO human annotations for go_id + descendants (TSV as DataFrame). Cached to
    cache_path if given (QuickGO is slow, and this pins the gene set for reproducibility)."""
    if cache_path is not None and os.path.exists(cache_path):
        return pd.read_csv(cache_path, sep="\t")

    import requests
    r = requests.get(QUICKGO_URL.format(go_id=go_id), headers={"Accept": "text/tsv"}, timeout=timeout)
    r.raise_for_status()
    ann = pd.read_csv(io.StringIO(r.text), sep="\t")
    if cache_path is not None:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        ann.to_csv(cache_path, sep="\t", index=False)
    return ann


def go_gene_symbols(ann):
    """Clean HGNC symbols from a QuickGO annotation table ("involved_in" only)."""
    ann = ann[ann["QUALIFIER"] == "involved_in"]
    symbols = ann["SYMBOL"].replace(SYMBOL_ALIASES)
    is_accession = symbols.str.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]")
    symbols = symbols[~is_accession & ~symbols.str.startswith("DKFZ")]
    return sorted(symbols.unique())


def gene_set_long(genes, df_genes, meta, compartment="Nuclear", pseudocount=1.0):
    """Long table (one row per gene x sample) for the genes present in df_genes, restricted
    to one compartment, with log2fc (vs. own D0), log2_expression and a per-gene z-score of
    log2_expression across all samples of that compartment."""
    present = [g for g in genes if g in df_genes.index]
    long = pd.concat(
        [gene_expression_over_time(g, df_genes, meta, pseudocount=pseudocount).assign(gene=g) for g in present],
        ignore_index=True,
    )
    long = long[long["Compartment"] == compartment].copy()
    grp = long.groupby("gene")["log2_expression"]
    long["zscore"] = (long["log2_expression"] - grp.transform("mean")) / grp.transform("std")
    return long


def pathway_score_per_sample(long):
    """One row per sample: median log2FC vs. D0 and mean z-score across the gene set."""
    return (
        long.groupby(["Sample", "Day", "Genotype"])
        .agg(median_log2fc=("log2fc", "median"), mean_zscore=("zscore", "mean"), n_genes=("gene", "nunique"))
        .reset_index()
    )


def _spearman_vs_day(df, value):
    rho, p = spearmanr(df["Day"], df[value])
    return pd.Series({"rho": rho, "p": p, "n": len(df)})


def pathway_trend(scores):
    """Spearman(score, Day) per genotype, for both pathway scores."""
    rows = []
    for value in ["median_log2fc", "mean_zscore"]:
        for genotype, sub in scores.groupby("Genotype"):
            rows.append({"score": value, "Genotype": genotype, **_spearman_vs_day(sub, value)})
    return pd.DataFrame(rows)


def gene_trends(long, value="log2_expression"):
    """Per gene x genotype: Spearman(value, Day), BH q within genotype, and the mean
    log2FC at each day (so the D35-vs-D0 amplitude sits next to the trend)."""
    trend = (
        long.groupby(["Genotype", "gene"])[["Day", value]]
        .apply(_spearman_vs_day, value=value)
        .reset_index()
    )
    trend["q"] = trend.groupby("Genotype")["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    trend["direction"] = np.select(
        [(trend["q"] < 0.05) & (trend["rho"] > 0), (trend["q"] < 0.05) & (trend["rho"] < 0)],
        ["up", "down"], default="n.s.",
    )
    mean_fc = (
        long.groupby(["Genotype", "gene", "Day"])["log2fc"].mean()
        .unstack("Day").add_prefix("mean_log2fc_D").reset_index()
    )
    return trend.merge(mean_fc, on=["Genotype", "gene"])


def plot_pathway_over_time(long, scores, metric="log2fc", title="", figsize=(9, 4.2)):
    """One panel per genotype (CTRL | VCP), shared y. Thin grey lines = per-gene mean across
    replicates; bold line = pathway score (median log2FC, or mean z-score), mean +- SD
    across replicates. metric: "log2fc" or "zscore"."""
    score_col = {"log2fc": "median_log2fc", "zscore": "mean_zscore"}[metric]
    ylabel = {"log2fc": "log2FC vs. D0", "zscore": "z-score of log2(expression + 1)"}[metric]
    score_label = {"log2fc": "median across genes", "zscore": "mean across genes"}[metric] + "\n(mean ± SD across replicates)"

    gene_means = long.groupby(["Genotype", "gene", "Day"])[metric].mean().reset_index()

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    for ax, genotype in zip(axes, GENOTYPE_ORDER):
        for _, g in gene_means[gene_means["Genotype"] == genotype].groupby("gene"):
            ax.plot(g["Day"], g[metric], color="lightgrey", linewidth=0.7, zorder=1)
        sns.lineplot(
            data=scores[scores["Genotype"] == genotype], x="Day", y=score_col,
            color=GENOTYPE_PALETTE[genotype], marker="o", linewidth=2.5, errorbar="sd", err_style="bars",
            err_kws={"capsize": 4, "elinewidth": 1.2}, label=score_label, ax=ax, zorder=3,
        )
        if metric == "log2fc":
            ax.axhline(0, color="grey", linestyle=":", linewidth=1)
        ax.set_xticks(DAY_ORDER)
        ax.set_xlabel("Day")
        ax.set_title(GENOTYPE_LABEL[genotype])
        ax.legend(frameon=False, loc="upper left")
        sns.despine(ax=ax)
    axes[0].set_ylabel(ylabel)
    axes[1].set_ylabel("")
    fig.suptitle(title, weight="bold")
    plt.tight_layout()
    return fig, axes


def plot_pathway_score_combined(scores, metric="log2fc", title="", figsize=(5.5, 4.2), ax=None):
    """Pathway score only (no per-gene lines), CTRL and VCP on the same axis, mean +- SD
    across replicates. metric: "log2fc" (median log2FC) or "zscore" (mean z-score)."""
    score_col = {"log2fc": "median_log2fc", "zscore": "mean_zscore"}[metric]
    ylabel = {
        "log2fc": "log2FC vs. D0\n(median across genes)",
        "zscore": "z-score of log2(expression + 1)\n(mean across genes)",
    }[metric]

    fig = ax.figure if ax is not None else plt.figure(figsize=figsize)
    ax = ax or plt.gca()
    plot_df = scores.assign(Genotype=scores["Genotype"].map(GENOTYPE_LABEL))
    sns.lineplot(
        data=plot_df, x="Day", y=score_col, hue="Genotype",
        hue_order=[GENOTYPE_LABEL[g] for g in GENOTYPE_ORDER],
        palette={GENOTYPE_LABEL[g]: c for g, c in GENOTYPE_PALETTE.items()},
        marker="o", linewidth=2.5, errorbar="sd", err_style="bars", err_kws={"capsize": 4, "elinewidth": 1.2}, ax=ax,
    )
    ax.axhline(0, color="grey", linestyle=":", linewidth=1)
    ax.set_xticks(DAY_ORDER)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.set_title(title, weight="bold")
    ax.legend(title="mean ± SD across replicates", title_fontsize=8, frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    return fig, ax


def plot_pathway_heatmap(long, metric="log2fc", title="", vlim=None, figsize=(7, 11),
                         row_order=None, cmap="PuOr_r", row_breaks=None, ytick_size=7, cbar_shrink=0.4):
    """Genes x Day heatmap, CTRL | VCP side by side, each cell = mean across replicates of
    `metric`: "log2fc" (vs. own D0) or "zscore" (per-gene z-score of log2(expression + 1)
    across all samples of the compartment, both genotypes). Rows are clustered once on the
    concatenated CTRL+VCP profiles, so both panels share the order -- unless `row_order`
    (list of genes) is given, then rows follow it as is. `row_breaks`: genes after which a
    white separator line is drawn (e.g. to split gene groups)."""
    cbar_label = {"log2fc": "mean log2FC vs. D0", "zscore": "mean z-score of log2(expression + 1)"}[metric]
    mat = long.groupby(["gene", "Genotype", "Day"])[metric].mean().unstack(["Genotype", "Day"])
    mat = mat[[(g, d) for g in GENOTYPE_ORDER for d in DAY_ORDER]]
    if row_order is None:
        row_order = mat.index[leaves_list(linkage(mat.values, method="average", metric="euclidean"))]
    mat = mat.loc[[g for g in row_order if g in mat.index]]
    vlim = vlim or np.nanpercentile(np.abs(mat.values), 98)

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    for i, (ax, genotype) in enumerate(zip(axes, GENOTYPE_ORDER)):
        sns.heatmap(
            mat[genotype], cmap=cmap, center=0, vmin=-vlim, vmax=vlim, ax=ax,
            cbar=i == 1, cbar_kws={"label": cbar_label, "shrink": cbar_shrink},
            yticklabels=True, linewidths=0,
        )
        for gene in row_breaks or []:
            ax.axhline(mat.index.get_loc(gene) + 1, color="white", linewidth=4)
        ax.set_title(GENOTYPE_LABEL[genotype])
        ax.set_xlabel("Day")
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=ytick_size, rotation=0)
    fig.suptitle(title, weight="bold")
    plt.tight_layout()
    return fig, axes, mat


def genotype_difference_per_day(long_gene, value="log2_expression"):
    """One gene, per Day: mean CTRL, mean VCP, difference VCP - CTRL with its Welch 95% CI,
    Welch t-test p (replicates of each genotype as independent samples) and BH q across days.
    With value="log2_expression" the difference is a log2 fold change VCP / CTRL."""
    rows = []
    for day, sub in long_gene.groupby("Day"):
        wt = sub.loc[sub["Genotype"] == "WT", value].values
        vcp = sub.loc[sub["Genotype"] == "VCP", value].values
        diff = vcp.mean() - wt.mean()
        v_wt, v_vcp = wt.var(ddof=1) / len(wt), vcp.var(ddof=1) / len(vcp)
        se = np.sqrt(v_wt + v_vcp)
        dof = (v_wt + v_vcp) ** 2 / (v_wt ** 2 / (len(wt) - 1) + v_vcp ** 2 / (len(vcp) - 1))
        half = t_dist.ppf(0.975, dof) * se
        rows.append({
            "Day": day, "mean_CTRL": wt.mean(), "mean_VCP": vcp.mean(), "diff_VCP_minus_CTRL": diff,
            "ci_low": diff - half, "ci_high": diff + half,
            "p_welch": ttest_ind(vcp, wt, equal_var=False).pvalue, "n_CTRL": len(wt), "n_VCP": len(vcp),
        })
    res = pd.DataFrame(rows)
    res["q_BH"] = multipletests(res["p_welch"], method="fdr_bh")[1]
    return res


def _stars(q):
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def plot_gene_genotype_difference(long_gene, gene, value="log2_expression", title="", figsize=(9.5, 4)):
    """Left: the gene over time, CTRL vs. VCP (mean +- SD across replicates, replicates as
    dots). Right: VCP - CTRL per day with Welch 95% CI; stars = BH q across days
    (* < 0.05, ** < 0.01, *** < 0.001). Returns (fig, axes, difference table)."""
    ylabel = {"log2_expression": "log2(expression + 1)", "zscore": "z-score of log2(expression + 1)",
              "log2fc": "log2FC vs. D0"}[value]
    res = genotype_difference_per_day(long_gene, value=value)
    plot_df = long_gene.assign(Genotype=long_gene["Genotype"].map(GENOTYPE_LABEL))
    palette = {GENOTYPE_LABEL[g]: c for g, c in GENOTYPE_PALETTE.items()}
    hue_order = [GENOTYPE_LABEL[g] for g in GENOTYPE_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    ax = axes[0]
    sns.lineplot(
        data=plot_df, x="Day", y=value, hue="Genotype", hue_order=hue_order, palette=palette,
        marker="o", linewidth=2.5, errorbar="sd", err_style="bars",
        err_kws={"capsize": 4, "elinewidth": 1.2}, ax=ax,
    )
    for genotype, sub in plot_df.groupby("Genotype"):
        offset = -0.6 if genotype == "CTRL" else 0.6
        ax.scatter(sub["Day"] + offset, sub[value], s=10, color=palette[genotype], alpha=0.5, linewidths=0)
    ax.set_xticks(DAY_ORDER)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.legend(title="mean ± SD across replicates", title_fontsize=8, frameon=False)
    ax.set_title("Expression over time")
    sns.despine(ax=ax)

    ax = axes[1]
    ax.errorbar(
        res["Day"], res["diff_VCP_minus_CTRL"],
        yerr=[res["diff_VCP_minus_CTRL"] - res["ci_low"], res["ci_high"] - res["diff_VCP_minus_CTRL"]],
        color=GENOTYPE_PALETTE["VCP"], marker="o", linewidth=2, capsize=4,
    )
    ax.axhline(0, color="grey", linestyle=":", linewidth=1)
    for _, r in res.iterrows():
        ax.text(r["Day"], r["ci_high"], _stars(r["q_BH"]), ha="center", va="bottom", fontsize=11)
    ax.set_xticks(DAY_ORDER)
    ax.set_xlabel("Day")
    ax.set_ylabel(f"VCP - CTRL\n({ylabel})")
    ax.set_title("VCP - CTRL (mean, Welch 95% CI)")
    sns.despine(ax=ax)

    fig.suptitle(title or gene, weight="bold")
    plt.tight_layout()
    return fig, axes, res


def safe_go_id(go_id):
    return re.sub(r"[^A-Za-z0-9]", "_", go_id)
