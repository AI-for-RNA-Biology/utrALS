"""Plot a gene's expression over the differentiation time course (D0-D35), either as log2
fold change vs. D0 or as plain log2(expression), colored by genotype with compartment shown
as line style (solid = Nuclear, dashed = Cytoplasmic) on a single axis.

Source: df_genes_filtered_normalized.csv, one column per sample named
"{Sample}_{Day}_{Compartment}_{Genotype}" (e.g. "CTRL1_D0_Nuclear_Wildtype",
"MUT2_D22_Cytoplasmic_VCP-R191Q"). log2FC is computed per (Sample, Compartment) trajectory
against that same replicate/compartment's own D0 value, so each biological replicate is its
own baseline -- not pooled/averaged before the log2 transform.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DAY_ORDER = [0, 3, 7, 14, 22, 35]
GENOTYPE_PALETTE = {"WT": "black", "VCP": "#9047d6"}
COMPARTMENT_ORDER = ["Nuclear", "Cytoplasmic"]
COMPARTMENT_DASHES = {"Nuclear": "", "Cytoplasmic": (4, 2)}


def load_expression_metadata(expr_path):
    """Load df_genes_filtered_normalized.csv and parse its sample columns into
    (Sample, Day, Compartment, Genotype). Returns (df_genes, metadata).

    The first column is used as the gene index ("geneName" here), so any table with the
    same "{Sample}_{Day}_{Compartment}_{Genotype}" sample columns also loads (e.g. the
    GSE152983 proteomics, indexed by "gene_name")."""
    df_genes = pd.read_csv(expr_path, index_col=0)

    sample_cols = df_genes.columns.tolist()
    meta = pd.DataFrame({"sample": sample_cols})
    parts = meta["sample"].str.split("_", n=3, expand=True)
    meta["Sample"] = parts[0]
    meta["Day"] = parts[1].str[1:].astype(int)
    meta["Compartment"] = parts[2]
    meta["Genotype"] = parts[3].replace({
        "Wildtype": "WT", "VCP-R155C": "VCP", "VCP-R191Q": "VCP",
    })
    meta = meta.set_index("sample")
    return df_genes, meta


def gene_expression_over_time(gene, df_genes, meta, pseudocount=1.0):
    """Long-format table for one gene: sample, Sample, Day, Compartment, Genotype,
    expression, log2_expression, log2fc (vs. that same Sample+Compartment's own D0 value)."""
    if gene not in df_genes.index:
        raise KeyError(f"{gene!r} not found in expression matrix")

    expr = df_genes.loc[gene]
    if isinstance(expr, pd.DataFrame):  # duplicate gene rows -- sum them
        expr = expr.sum(axis=0)

    long = meta.copy()
    long["expression"] = expr.reindex(long.index)
    long["log2_expression"] = np.log2(long["expression"] + pseudocount)

    baseline = (
        long[long["Day"] == 0]
        .set_index(["Sample", "Compartment"])["expression"]
    )
    long["baseline"] = long.set_index(["Sample", "Compartment"]).index.map(baseline)
    long["log2fc"] = np.log2(long["expression"] + pseudocount) - np.log2(long["baseline"] + pseudocount)
    return long.reset_index(drop=True)


# kept as an alias -- this function used to be the only one, now a thin wrapper
gene_log2fc_over_time = gene_expression_over_time

_METRIC_YLABELS = {
    "log2fc": "log2FC vs. D0 (pseudocount={pc:g})",
    "log2_expression": "log2(expression + {pc:g})",
}


def plot_gene_over_time(gene, df_genes, meta, metric="log2fc", pseudocount=1.0, figsize=(6, 4.5), ax=None):
    """Single-axis line plot over the time course: color = genotype, line style = compartment
    (solid = Nuclear, dashed = Cytoplasmic). `metric` is "log2fc" (vs. that replicate's own D0,
    default) or "log2_expression" (plain log2 of the normalized counts, no D0 baseline).
    Error bars are 95% CI across replicates within each genotype x compartment x day.
    Returns (fig, ax, long_df)."""
    if metric not in _METRIC_YLABELS:
        raise ValueError(f"metric must be one of {list(_METRIC_YLABELS)}, got {metric!r}")

    long = gene_expression_over_time(gene, df_genes, meta, pseudocount=pseudocount)

    fig = ax.figure if ax is not None else plt.figure(figsize=figsize)
    ax = ax or plt.gca()

    sns.lineplot(
        data=long, x="Day", y=metric, hue="Genotype", hue_order=["WT", "VCP"],
        style="Compartment", style_order=COMPARTMENT_ORDER, dashes=COMPARTMENT_DASHES,
        palette=GENOTYPE_PALETTE, marker="o", err_style="bars",
        err_kws={"capsize": 4, "elinewidth": 1}, ax=ax,
    )
    if metric == "log2fc":
        ax.axhline(0, color="grey", linestyle=":", linewidth=1)
    ax.set_xticks(DAY_ORDER)
    ax.set_xlabel("Day")
    ax.set_ylabel(_METRIC_YLABELS[metric].format(pc=pseudocount))
    ax.set_title(gene, weight="bold")
    sns.despine(ax=ax)
    plt.tight_layout()
    return fig, ax, long


def plot_gene_log2fc_over_time(gene, df_genes, meta, pseudocount=1.0, figsize=(9, 4.5)):
    """Two-panel (Nuclear | Cytoplasmic) version, kept for compatibility with the earlier
    faceted layout -- prefer plot_gene_over_time() for the combined single-axis view."""
    long = gene_expression_over_time(gene, df_genes, meta, pseudocount=pseudocount)

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    for ax, compartment in zip(axes, COMPARTMENT_ORDER):
        sub = long[long["Compartment"] == compartment]
        sns.lineplot(
            data=sub, x="Day", y="log2fc", hue="Genotype", hue_order=["WT", "VCP"],
            palette=GENOTYPE_PALETTE, marker="o", err_style="bars",
            err_kws={"capsize": 4, "elinewidth": 1}, ax=ax,
        )
        ax.axhline(0, color="grey", linestyle=":", linewidth=1)
        ax.set_xticks(DAY_ORDER)
        ax.set_xlabel("Day")
        ax.set_title(compartment)
        sns.despine(ax=ax)

    axes[0].set_ylabel(f"log2FC vs. D0 (pseudocount={pseudocount:g})")
    axes[1].set_ylabel("")
    fig.suptitle(gene, weight="bold")
    plt.tight_layout()
    return fig, axes, long


# ── Gene / protein expression over time (WT vs VCP heatmaps) ──────────────────
import matplotlib.patches as mpatches
from scipy.stats import ttest_ind
from scipy.cluster.hierarchy import linkage, leaves_list

from .rbp_names import RBP_TO_SYMBOL


def expression_long(rbp_list, df, meta, value_col, compartment="Nuclear"):
    """Per-replicate long table for each RBP (one compartment), built with
    gene_expression_over_time; rows labelled with the CLIP RBP name. Exact symbol match.
    `df, meta` come from gene_expression_over_time.load_expression_metadata."""
    frames, missing = [], []
    for rbp in dict.fromkeys(rbp_list):   # drop duplicate names, keep order
        symbol = RBP_TO_SYMBOL.get(rbp, rbp)
        if symbol not in df.index:
            missing.append(rbp)
            continue
        long = gene_expression_over_time(symbol, df, meta)
        long = long[long["Compartment"] == compartment]
        frames.append(long.assign(RBP=rbp, value=long[value_col]))
    if missing:
        print(f"  not found: {missing}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def plot_wt_vs_vcp_heatmap(rbp_list, title, save_path, df, meta, value_col="log2_expression",
                           compartment="Nuclear"):
    """Expression heatmap over time, WT vs VCP side by side. Each RBP row is z-scored jointly
    over the 12 columns (6 days x WT + VCP), so colours are comparable between the two
    panels (level) as well as across days (trend). Welch's t-test (WT vs VCP, per RBP per day, on the
    per-replicate values) marking significant days with '*'.

    value_col: "log2_expression" for RNA (log2(normalised counts + 1)), "expression" for the
    proteomics (already log2), or "log2fc" (vs. each replicate's own D0)."""
    long = expression_long(rbp_list, df, meta, value_col, compartment)
    if long.empty:
        return
    means = (long.groupby(["Genotype", "RBP", "Day"])["value"].mean()
             .unstack("Day").reindex(columns=DAY_ORDER))
    wt, vcp = means.loc["WT"], means.loc["VCP"]
    common = wt.index.intersection(vcp.index)
    common = common[(wt.loc[common].notna().sum(axis=1) >= 2) & (vcp.loc[common].notna().sum(axis=1) >= 2)]
    no_data = [rbp for rbp in dict.fromkeys(long["RBP"]) if rbp not in common]
    if no_data:
        print(f"  no {compartment.lower()} values over time: {no_data}")
    wt, vcp = wt.loc[common], vcp.loc[common]
    if wt.empty:
        return

    both = pd.concat([wt, vcp], axis=1)
    row_mean, row_std = both.mean(axis=1), both.std(axis=1)
    wt_scaled = wt.sub(row_mean, axis=0).div(row_std, axis=0)
    vcp_scaled = vcp.sub(row_mean, axis=0).div(row_std, axis=0)

    order = (wt_scaled.index[leaves_list(linkage(wt_scaled.fillna(0), method="average"))]
             if len(wt_scaled) > 1 else wt_scaled.index)
    wt_scaled, vcp_scaled = wt_scaled.loc[order], vcp_scaled.loc[order]

    sig = pd.DataFrame(False, index=order, columns=DAY_ORDER)
    for (rbp, day), grp in long[long["RBP"].isin(order)].groupby(["RBP", "Day"]):
        wt_vals = grp.loc[grp["Genotype"] == "WT", "value"]
        vcp_vals = grp.loc[grp["Genotype"] == "VCP", "value"]
        _, p = ttest_ind(wt_vals, vcp_vals, equal_var=False, nan_policy="omit")
        sig.loc[rbp, day] = p < 0.05

    vmax = np.nanmax(np.abs(np.concatenate([wt_scaled.values, vcp_scaled.values])))
    # Colorbar in its own narrow column, so the WT and VCP panels get the same width
    fig, (ax_wt, ax_vcp, cax) = plt.subplots(1, 3, figsize=(10.5, 0.4 * len(order) + 2),
                                             gridspec_kw={"width_ratios": [1, 1, 0.05]})
    axes = [ax_wt, ax_vcp]
    sns.heatmap(wt_scaled, cmap="PuOr_r", center=0, vmin=-vmax, vmax=vmax, ax=ax_wt, cbar=False)
    sns.heatmap(vcp_scaled, cmap="PuOr_r", center=0, vmin=-vmax, vmax=vmax, ax=ax_vcp,
                cbar_ax=cax, cbar_kws={"label": "z-score (row, WT + VCP)"})
    ax_vcp.set_yticks([])   # same rows as WT: labels only on the left panel
    for ax, genotype in zip(axes, ["WT", "VCP"]):
        ax.set_facecolor("lightgrey")   # NaN (e.g. protein not detected)
        ax.set_title(genotype)
        ax.set_xticklabels([f"D{d}" for d in DAY_ORDER], rotation=0)
        ax.set_xlabel("Day")
        ax.set_ylabel("")
        for i, rbp in enumerate(order):
            for j, day in enumerate(DAY_ORDER):
                if sig.loc[rbp, day]:
                    ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", color="black")

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight", metadata={"CreationDate": None} if str(save_path).endswith(".pdf") else None)  # no date: identical file on re-runs
    plt.show()


def plot_expression_at_day(region_rbp_lists, df, meta, value_col="log2_expression", compartment="Nuclear",
                           day=35, value_label="log2(expression + 1)", title=None, save_path=None,
                           palette=None):
    """Compare expression LEVELS across RBPs at one day (no row normalisation).

    One horizontal panel per region list (e.g. {"R1": [...], "R2": [...], "R3": [...]}),
    sharing the x-axis. Each RBP: mean WT and VCP bars + per-replicate dots, RBPs sorted by
    WT mean (most expressed on top). '*' after the RBP name: Welch t-test WT vs VCP p < 0.05.
    Returns the long per-replicate table used for the plot."""
    palette = palette or {"WT": "#c5c5c5", "VCP": "#9047d6"}
    longs = {}
    for region, rbps in region_rbp_lists.items():
        long = expression_long(rbps, df, meta, value_col, compartment)
        longs[region] = long[long["Day"] == day].dropna(subset=["value"]) if not long.empty else long
    longs = {r: l for r, l in longs.items() if not l.empty}
    if not longs:
        return pd.DataFrame()

    n_max = max(l["RBP"].nunique() for l in longs.values())
    fig, axes = plt.subplots(1, len(longs), figsize=(4.2 * len(longs), 0.35 * n_max + 1.8),
                             sharex=True, squeeze=False)
    for ax, (region, long) in zip(axes[0], longs.items()):
        order = (long[long["Genotype"] == "WT"].groupby("RBP")["value"].mean()
                 .reindex(long["RBP"].unique()).sort_values(ascending=False).index.tolist())
        labels = []
        for rbp in order:
            sub = long[long["RBP"] == rbp]
            wt, vcp = sub.loc[sub["Genotype"] == "WT", "value"], sub.loc[sub["Genotype"] == "VCP", "value"]
            p = ttest_ind(wt, vcp, equal_var=False).pvalue if len(wt) > 1 and len(vcp) > 1 else np.nan
            labels.append(f"{rbp} *" if p < 0.05 else rbp)
        # ci=None / no legend= kwarg: compatible with seaborn 0.11 (utrals environment)
        sns.barplot(data=long, y="RBP", x="value", hue="Genotype", order=order, hue_order=["WT", "VCP"],
                    palette=palette, ci=None, orient="h", ax=ax)
        sns.stripplot(data=long, y="RBP", x="value", hue="Genotype", order=order, hue_order=["WT", "VCP"],
                      dodge=True, palette={"WT": "black", "VCP": "black"}, size=2.5, orient="h", ax=ax)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(region, weight="bold")
        ax.set_xlabel(value_label)
        ax.set_ylabel("")
        sns.despine(ax=ax)
    for ax in axes[0][len(longs):]:
        ax.axis("off")
    axes[0][len(longs) - 1].legend(handles=[mpatches.Patch(color=palette[g], label=g) for g in ["WT", "VCP"]],
                                   frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(title or f"Expression at D{day} ({compartment})", weight="bold")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", metadata={"CreationDate": None} if str(save_path).endswith(".pdf") else None)  # no date: identical file on re-runs
    plt.show()
    return pd.concat(longs, names=["region"]).reset_index(level=0)
