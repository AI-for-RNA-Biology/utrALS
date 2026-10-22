"""Welch profiles across the MN differentiation time course (MN_diff): per day, genotype and fraction, at the
proximal and distal PAS and around the stop codon; by-region zooms; Fisher profiles of consecutive days; Fisher heatmaps.

Layout `time_course` of scripts/clip/clip_analysis.py.
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from .panels import RC, LINEWIDTH, panel_grid, block_title, figure_title
from .fisher import compute_signed_enrichment, time_course_title
from .common import setup_output_directories, load_integrated, save, for_each, fisher_heatmaps


palette_days = {
    "D0": "#1f77b4",
    "D3": "#0ea3ff",
    "D7": "#c4dbdf",
    "D14": "#a5baa3c5",
    "D22": "#67b46b",
    "D35": "#126e02"
}


palette_mutations = {
    "Wildtype": "#2d609a",
    "VCP": "#d17fb0"
}


def _draw_welch_profiles_across_time(axes, df_results, RBP, treatments, palette_days,
                                     anchor='pPAS', sign=1, comparison_day='D35',
                                     nt_before=-1000, length_region=50):
    """
    Draw Welch profiles (lines only) onto a 3x2 axes grid (columns: Nuclear, Cytoplasmic).
    Rows 1-2: all days for Wildtype / VCP. Row 3: `comparison_day` Wildtype vs VCP.
    `sign` flips the profile (distal profiles are plotted as -value).
    `nt_before` / `length_region`: window used in the integration (x = start of each region).
    Labels are kept on the outer panels only; legends are left to the caller.
    """
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in treatments]
    # pvalue == 0 underflows -log10(pvalue) to +-inf; cap it at the most extreme finite
    # value actually observed in this RBP's data instead of a fixed constant, so a single
    # degenerate (tiny-sample) bin doesn't blow up the shared y-axis.
    raw = df_results[RBP][cols].apply(pd.to_numeric, errors='coerce')
    finite = raw.replace([np.inf, -np.inf], np.nan)
    clipped = raw.replace([np.inf, -np.inf], [finite.max().max(), finite.min().min()]) * sign

    ylim = (clipped.min().min() - 2, clipped.max().max() + 2)
    distances = [nt_before + i * length_region for i in range(df_results[RBP].shape[0])]

    def style(ax, title, row, col):
        ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
        if row == axes.shape[0] - 1:
            ax.set_xlabel(f'Distance to {anchor} (nt)')
        if col == 0:
            ax.set_ylabel('-log10(pval)*dPUD')
        ax.set_title(title)
        ax.set_ylim(ylim)
        sns.despine(ax=ax)

    for row, mutation in enumerate(["Wildtype", "VCP"]):
        for col, compartment in enumerate(["Nuclear", "Cytoplasmic"]):
            ax = axes[row, col]
            for treatment in treatments:
                if compartment in treatment and mutation in treatment:
                    day = treatment.split('_')[0]
                    ax.plot(distances, clipped[f'-log10(pval)*dPUD_{treatment}'], label=day,
                            color=palette_days[day], linewidth=LINEWIDTH)
            style(ax, f'{compartment} {mutation}', row, col)

    for col, compartment in enumerate(["Nuclear", "Cytoplasmic"]):
        ax = axes[2, col]
        for mutation in ["Wildtype", "VCP"]:
            ax.plot(distances, clipped[f'-log10(pval)*dPUD_{comparison_day}_{compartment}_{mutation}'],
                    label=mutation, color=palette_mutations[mutation], linewidth=LINEWIDTH)
        style(ax, f'{compartment} {comparison_day} Wildtype vs VCP', 2, col)


def _welch_profiles_figure(blocks, RBP, treatments, palette_days, comparison_day='D35'):
    """
    One 3x2 block of PANEL_CM x PANEL_CM Welch profile panels per entry of `blocks`
    (df_results, anchor, sign, name, nt_before, length_region), side by side, with the
    legends at the bottom. Save with bbox_inches='tight'.
    """
    with plt.rc_context(RC):
        fig, axes = panel_grid(3, 2 * len(blocks), block_size=2, top_cm=1.6, bottom_cm=1.6)
        for b, (df_results, anchor, sign, name, nt_before, length_region) in enumerate(blocks):
            block_axes = axes[:, 2 * b:2 * b + 2]
            _draw_welch_profiles_across_time(block_axes, df_results, RBP, treatments, palette_days,
                                             anchor=anchor, sign=sign, comparison_day=comparison_day,
                                             nt_before=nt_before, length_region=length_region)
            block_title(fig, block_axes[0, 0], block_axes[0, 1], name)
        figure_title(fig, axes, RBP)

        day_handles, day_labels = axes[0, 0].get_legend_handles_labels()
        mut_handles, mut_labels = axes[2, 0].get_legend_handles_labels()
        blank = plt.Line2D([], [], linestyle='none')
        fig.legend([blank] + day_handles + [blank] + mut_handles,
                   ['Days:'] + day_labels + [f'   {comparison_day}:'] + mut_labels,
                   loc='lower center', ncol=len(day_labels) + len(mut_labels) + 2,
                   frameon=False, handlelength=1.5, columnspacing=1.2)
    return fig


def welch_profiles_across_time_along_prox(df_results_proximals, RBP, treatments, palette_days,
                                          nt_before=-1000, length_region=50):
    return _welch_profiles_figure([(df_results_proximals, 'pPAS', 1, 'Proximal', nt_before, length_region)],
                                  RBP, treatments, palette_days)


def welch_profiles_across_time_along_dist(df_results_distals, RBP, treatments, palette_days,
                                          nt_before=-1000, length_region=50):
    return _welch_profiles_figure([(df_results_distals, 'dPAS', -1, 'Distal', nt_before, length_region)],
                                  RBP, treatments, palette_days)


def welch_profiles_across_time_combined(df_results_proximals, df_results_distals, RBP, treatments, palette_days,
                                        df_results_stop_codon=None, nt_before=-1000, length_region=50,
                                        stop_codon_nt_before=-100, stop_codon_length_region=25):
    """
    Stop codon (left, if `df_results_stop_codon` has this RBP), proximal and distal (right)
    3x2 Welch profile blocks on a single page, with shared legends at the bottom.
    Stop codon profiles keep their sign, like the proximal ones (positive = associated with
    proximal usage). Save with bbox_inches='tight'.
    """
    blocks = [(df_results_proximals, 'pPAS', 1, 'Proximal', nt_before, length_region),
              (df_results_distals, 'dPAS', -1, 'Distal', nt_before, length_region)]
    if df_results_stop_codon is not None and RBP in df_results_stop_codon:
        blocks.insert(0, (df_results_stop_codon, 'stop codon', 1, 'Stop codon',
                          stop_codon_nt_before, stop_codon_length_region))
    return _welch_profiles_figure(blocks, RBP, treatments, palette_days)


def welch_profiles_by_region(df_results_proximals, df_results_distals, RBP, treatments, palette_days,
                              positive_proximal_region=(-350, 0),
                              negative_proximal_region=(0, 200),
                              positive_distal_region=(-350, 0)):
    """
    Plot Welch profiles zoomed into specific regions around the proximal and distal PAS.

    Produces a 4 (Nuclear/Cytoplasmic x Wildtype/VCP) x 3 (region) grid of subplots.
    Each region column shares its own y-scale (computed across the 4 conditions for
    that region only), since the three regions are not expected to have comparable
    magnitudes.

    Parameters
    ----------
    df_results_proximals, df_results_distals : dict
        Dictionaries of DataFrames with results, keyed by RBP.
    RBP : str
    treatments : list
    palette_days : dict
    positive_proximal_region, negative_proximal_region, positive_distal_region : tuple
        (start, end) in nt around the corresponding PAS, defining each region's x-range.
    """
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in treatments]

    def clip_and_distances(df_results, region_size):
        raw = df_results[RBP][cols].apply(pd.to_numeric, errors='coerce')
        finite = raw.replace([np.inf, -np.inf], np.nan)
        clipped = raw.replace([np.inf, -np.inf], [finite.max().max(), finite.min().min()])
        distances = np.array([-1000 + i * 50 for i in range(df_results[RBP].shape[0])])
        return clipped, distances

    clipped_prox, distances_prox = clip_and_distances(df_results_proximals, df_results_proximals[RBP].shape[0])
    clipped_dist, distances_dist = clip_and_distances(df_results_distals, df_results_distals[RBP].shape[0])

    regions = [
        {'name': 'positive_proximal_region', 'clipped': clipped_prox, 'distances': distances_prox,
         'xrange': positive_proximal_region, 'xlabel': 'Distance to pPAS (nt)', 'sign': 1},
        {'name': 'negative_proximal_region', 'clipped': clipped_prox, 'distances': distances_prox,
         'xrange': negative_proximal_region, 'xlabel': 'Distance to pPAS (nt)', 'sign': 1},
        {'name': 'positive_distal_region', 'clipped': clipped_dist, 'distances': distances_dist,
         'xrange': positive_distal_region, 'xlabel': 'Distance to dPAS (nt)', 'sign': -1},
    ]

    conditions = [(mutation, compartment) for mutation in ["Wildtype", "VCP"] for compartment in ["Nuclear", "Cytoplasmic"]]

    # For each region, gather the per-condition day-values restricted to that region's
    # x-range, and the shared y-limits across all 4 conditions in that region.
    region_data = []
    for region in regions:
        xmin, xmax = region['xrange']
        mask = (region['distances'] >= xmin) & (region['distances'] <= xmax)
        x_masked = region['distances'][mask]

        per_condition = []
        region_values = []
        for mutation, compartment in conditions:
            sub_df = region['clipped'][[f'-log10(pval)*dPUD_{treatment}' for treatment in treatments
                                         if compartment in treatment and mutation in treatment]].copy()
            sub_df.columns = [treatment.split('_')[1] for treatment in sub_df.columns]
            sub_df = sub_df.iloc[mask.tolist()] * region['sign']
            per_condition.append(sub_df)
            region_values.append(sub_df.values)

        all_values = np.concatenate([v.ravel() for v in region_values])
        all_values = all_values[np.isfinite(all_values)]
        ymin, ymax = all_values.min() - 2, all_values.max() + 2
        region_data.append({'x': x_masked, 'per_condition': per_condition, 'ylim': (ymin, ymax)})

    with plt.rc_context(RC):
        # two-line panel titles: more space between rows
        fig, axes = panel_grid(len(conditions), len(regions), hspace_cm=1.3, top_cm=1.6)
        for row, (mutation, compartment) in enumerate(conditions):
            for col, region in enumerate(regions):
                ax = axes[row, col]
                sub_df = region_data[col]['per_condition'][row]
                x = region_data[col]['x']
                for day in sub_df.columns:
                    ax.plot(x, sub_df[day], label=day, color=palette_days[day], linewidth=LINEWIDTH)
                ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
                ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
                if row == len(conditions) - 1:
                    ax.set_xlabel(regions[col]['xlabel'])
                if col == 0:
                    ax.set_ylabel('-log10(pval)*dPUD')
                ax.set_title(f'{compartment} {mutation}\n{regions[col]["name"]}')
                ax.set_ylim(region_data[col]['ylim'])
                sns.despine(ax=ax)
                if row == 0 and col == len(regions) - 1:
                    ax.legend(title='Days', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

        figure_title(fig, axes, RBP, offset_cm=0.9)
    return fig


def plot_fisher_profiles_consecutive_days(RBP, 
                                          df_results, 
                                          days_pairs, 
                                          palette_days, 
                                          shift_type='distal',
                                          save_path=None):
    """
    Plot Fisher enrichment test profiles for consecutive day pairs across conditions.
    
    Parameters:
    -----------
    RBP : str
        RBP name to plot
    df_results_proximals : dict
        Dictionary of DataFrames containing fisher test results
    days_pairs : list
        List of tuples with consecutive day pairs (e.g., [("D0", "D3"), ("D3", "D7")])
    palette_days : dict
        Color palette for days
    save_path : str, optional
        Path to save the figure
    """
    # Calculate global min/max across all conditions
    all_enrichments = []
    for mutation in ["Wildtype", "VCP"]:
        for compartment in ["Nuclear", "Cytoplasmic"]:
            for pair in days_pairs:
                pairs = f"pairs_{shift_type}s_{pair[0]}_{compartment}_{mutation}_vs_{pair[1]}_{compartment}_{mutation}"
                signed_enrichment = compute_signed_enrichment(
                    df_results[RBP][f"fisher_results_{pairs}"]
                )
                all_enrichments.extend(signed_enrichment.values)
    
    global_ymin = min(all_enrichments) - 2
    global_ymax = max(all_enrichments) + 2
    
    # Create plots
    plt.figure(figsize=(10, 6))
    i = 1
    for mutation in ["Wildtype", "VCP"]:
        for compartment in ["Nuclear", "Cytoplasmic"]:
            plt.subplot(2, 2, i)
            
            for pair in days_pairs:
                pairs = f"pairs_{shift_type}s_{pair[0]}_{compartment}_{mutation}_vs_{pair[1]}_{compartment}_{mutation}"
                
                signed_enrichment = compute_signed_enrichment(
                    df_results[RBP][f"fisher_results_{pairs}"]
                )
                
                distances = [-1000+j*50 for j in range(df_results[RBP].shape[0])]
                plt.plot(distances, signed_enrichment, label=f"{pair[1]} vs {pair[0]}", 
                        color=palette_days[pair[1]], marker='o', markersize=3)
            
            plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
            plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
            plt.xlabel('Distance to pPAS (nt)')
            plt.ylabel("-log10(p-value)*sign")
            plt.title(f'{compartment} {mutation}')
            plt.ylim(global_ymin, global_ymax)
            
            if i == 2:
                plt.legend(title='Day pairs', bbox_to_anchor=(1.05, 1.2), loc='upper left')
            else:
                plt.legend().remove()
            
            sns.despine()
            i += 1
    
    plt.suptitle(f"Fisher enrichment test for {shift_type} pairs - {RBP}", weight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", metadata={"CreationDate": None} if str(save_path).endswith(".pdf") else None)  # no date: identical file on re-runs



def run(config, welch_only=False):
    """All outputs of one analysis config (see pipeline/datasets/MN_diff/configs/*analysis*.yaml). `welch_only`: only
    the Welch profiles across time (proximal, distal, combined); skip the by-region profiles and all Fisher outputs."""
    out = setup_output_directories(config['output_directory'], {
        'welch_proximal': 'welch_profiles_proximal', 'welch_distal': 'welch_profiles_distal',
        'welch_combined': 'welch_profiles_combined', 'welch_regions': 'welch_profiles_by_region',
        'fisher_consecutive_proximal': 'fisher_consecutive_days_proximal',
        'fisher_consecutive_distal': 'fisher_consecutive_days_distal',
        'fisher_heatmaps_proximal': 'fisher_heatmaps_proximal', 'fisher_heatmaps_distal': 'fisher_heatmaps_distal'})
    proximals, distals, stop_codon = load_integrated(config)
    rbps_proximal, rbps_distal = list(proximals), list(distals)
    rbps_common = [rbp for rbp in rbps_proximal if rbp in distals]
    treatments, palette = config['treatments'], config['palette_days']

    for_each(rbps_proximal, "Welch profiles (proximal)", lambda rbp: (
        welch_profiles_across_time_along_prox(proximals, rbp, treatments, palette),
        save(out['welch_proximal'] / f"welch_profile_proximal_{rbp}.pdf")))
    for_each(rbps_distal, "Welch profiles (distal)", lambda rbp: (
        welch_profiles_across_time_along_dist(distals, rbp, treatments, palette),
        save(out['welch_distal'] / f"welch_profile_distal_{rbp}.pdf")))
    # stop codon (when integrated for this RBP) + proximal + distal on one page
    window = dict(nt_before=config.get('nb_nucleotides_upstream', -1000), length_region=config.get('length_region', 50),
                  stop_codon_nt_before=config.get('stop_codon_nb_nucleotides_upstream', -100),
                  stop_codon_length_region=config.get('stop_codon_length_region', 25))
    for_each(rbps_common, "Welch profiles (stop codon + proximal + distal combined)", lambda rbp: (
        welch_profiles_across_time_combined(proximals, distals, rbp, treatments, palette,
                                            df_results_stop_codon=stop_codon, **window),
        save(out['welch_combined'] / f"{rbp}.pdf")))
    if welch_only:
        print("\n--welch-only: skipping by-region Welch profiles and Fisher outputs.")
        return

    for_each(rbps_common, "Welch profiles by region", lambda rbp: (
        welch_profiles_by_region(proximals, distals, rbp, treatments, palette),
        save(out['welch_regions'] / f"welch_profile_by_region_{rbp}.pdf")))
    days_pairs = [tuple(pair) for pair in config['days_pairs']]
    for iso, df_results, rbps in [('proximal', proximals, rbps_proximal), ('distal', distals, rbps_distal)]:
        for shift_type in ('distal', 'proximal'):
            for_each(rbps, f"Fisher consecutive days profiles ({iso}, {shift_type} shifts)", lambda rbp: (
                plot_fisher_profiles_consecutive_days(rbp, df_results, days_pairs, palette, shift_type=shift_type),
                save(out[f'fisher_consecutive_{iso}'] / f"fisher_consecutive_{iso}_{shift_type}_shift_{rbp}_.pdf")))
    fisher_heatmaps(config, {'proximal': proximals, 'distal': distals},
                    {'proximal': out['fisher_heatmaps_proximal'], 'distal': out['fisher_heatmaps_distal']},
                    title=time_course_title, merge_ptb=True)
