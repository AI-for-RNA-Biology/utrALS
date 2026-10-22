"""Welch profiles of a knockdown experiment (SFPQ_TDP_FUS_KD, one analysis per knockdown): rows = fractions, columns =
stop codon / proximal / distal, one line per treatment; Fisher heatmaps and Fisher profiles of the chosen pairs.

Layout `knockdown` of scripts/clip/clip_analysis.py.
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from .panels import RC, LINEWIDTH, panel_grid, figure_title
from .fisher import compute_signed_enrichment, treatment_title
from .common import setup_output_directories, load_integrated, save, for_each, fisher_heatmaps


def welch_profiles_2x2(df_results_proximals, df_results_distals, RBP, treatments, palette,
                       compartments={'Nuclear': '_Nuc', 'Cytoplasmic': '_Cyto'},
                       nb_nucleotides_upstream=-1000, length_region=50,
                       df_results_stop_codon=None, stop_codon_nb_nucleotides_upstream=-100,
                       stop_codon_length_region=25):
    """
    Welch profiles (lines only, panels.PANEL_CM panels) on a grid: rows = compartments
    (treatments ending with the given suffix), columns = stop codon (if `df_results_stop_codon`
    has this RBP) / proximal / distal. Distal profiles are plotted as -value; stop codon profiles
    keep their sign, like the proximal ones (positive = associated with proximal usage).
    Save with bbox_inches='tight'.
    """
    columns = [(df_results_proximals[RBP], 'pPAS', 1, 'Proximal', nb_nucleotides_upstream, length_region),
               (df_results_distals[RBP], 'dPAS', -1, 'Distal', nb_nucleotides_upstream, length_region)]
    if df_results_stop_codon is not None and RBP in df_results_stop_codon:
        columns.insert(0, (df_results_stop_codon[RBP], 'stop codon', 1, 'Stop codon',
                           stop_codon_nb_nucleotides_upstream, stop_codon_length_region))
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in treatments]

    with plt.rc_context(RC):
        fig, axes = panel_grid(len(compartments), len(columns), right_cm=4.0)
        for col, (df, anchor, sign, name, nt_before, region_length) in enumerate(columns):
            # pvalue == 0 gives +-inf: cap it at the most extreme finite value observed
            raw = df[cols].apply(pd.to_numeric, errors='coerce')
            finite = raw.replace([np.inf, -np.inf], np.nan)
            clipped = raw.replace([np.inf, -np.inf], [finite.max().max(), finite.min().min()]) * sign
            ylim = (clipped.min().min() - 2, clipped.max().max() + 2)
            distances = [nt_before + i * region_length for i in range(df.shape[0])]

            for row, (compartment, suffix) in enumerate(compartments.items()):
                ax = axes[row, col]
                for treatment in treatments:
                    if treatment.endswith(suffix):
                        ax.plot(distances, clipped[f'-log10(pval)*dPUD_{treatment}'], label=treatment,
                                color=palette[treatment], linewidth=LINEWIDTH)
                ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
                ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
                if row == len(compartments) - 1:
                    ax.set_xlabel(f'Distance to {anchor} (nt)')
                if col == 0:
                    ax.set_ylabel('-log10(pval)*dPUD')
                ax.set_title(f'{compartment} - {name}')
                ax.set_ylim(ylim)
                sns.despine(ax=ax)
                if col == len(columns) - 1:
                    ax.legend(title='Treatments', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

        figure_title(fig, axes, RBP, offset_cm=0.5)
    return fig


def plot_fisher_profiles(RBP, 
                         pairs_name,
                         df_results, 
                         palette_pairs, 
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
    for pair in pairs_name:
        signed_enrichment = compute_signed_enrichment(
            df_results[RBP][f"fisher_results_{pair}"]
        )
        all_enrichments.extend(signed_enrichment.values)
    
    global_ymin = min(all_enrichments) - 2
    global_ymax = max(all_enrichments) + 2
    
    # Create plots
    plt.figure(figsize=(12, 4))
        
    for pair in pairs_name:
        
        signed_enrichment = compute_signed_enrichment(
            df_results[RBP][f"fisher_results_{pair}"]
        )
        
        distances = [-1000+j*50 for j in range(df_results[RBP].shape[0])]
        plt.plot(distances, signed_enrichment, label=f"{pair}", 
                color=palette_pairs[pair], marker='o', markersize=3)
    
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.xlabel('Distance to pPAS (nt)')
    plt.ylabel("-log10(p-value)*sign")
    plt.title(f'{RBP}')
    plt.ylim(global_ymin, global_ymax)
    
    plt.legend(title='Pairs', bbox_to_anchor=(1.05, 1.2), loc='upper left')
    
    sns.despine()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')



def run(config, welch_only=False):
    """All outputs of one analysis config (see pipeline/datasets/SFPQ_TDP_FUS_KD/configs/*analysis*.yaml)."""
    out = setup_output_directories(config['output_directory'], {
        'welch': 'welch_profiles', 'fisher_heatmaps_proximal': 'fisher_heatmaps_proximal',
        'fisher_heatmaps_distal': 'fisher_heatmaps_distal', 'fisher_profiles_proximal': 'fisher_profiles_proximal',
        'fisher_profiles_distal': 'fisher_profiles_distal'})
    proximals, distals, stop_codon = load_integrated(config)
    rbps_common = [rbp for rbp in proximals if rbp in distals]
    treatments, palette = config['treatments'], config['palette']

    for_each(rbps_common, "Welch profiles (stop codon + proximal + distal)", lambda rbp: (
        welch_profiles_2x2(proximals, distals, rbp, treatments, palette,
                           nb_nucleotides_upstream=config.get('nb_nucleotides_upstream', -1000),
                           length_region=config.get('length_region', 50),
                           df_results_stop_codon=stop_codon,
                           stop_codon_nb_nucleotides_upstream=config.get('stop_codon_nb_nucleotides_upstream', -100),
                           stop_codon_length_region=config.get('stop_codon_length_region', 25)),
        save(out['welch'] / f"welch_profile_{rbp}.pdf")))
    if welch_only:
        return
    fisher_heatmaps(config, {'proximal': proximals, 'distal': distals},
                    {'proximal': out['fisher_heatmaps_proximal'], 'distal': out['fisher_heatmaps_distal']},
                    title=treatment_title, merge_ptb=False)
    pairs = config.get('pairs_for_profiles', [])
    if pairs:
        for iso, df_results in [('proximal', proximals), ('distal', distals)]:
            for_each(list(df_results), f"Fisher profiles ({iso})", lambda rbp: (
                plot_fisher_profiles(rbp, pairs, df_results, config.get('palette_pairs', {})),
                save(out[f'fisher_profiles_{iso}'] / f"fisher_profile_{iso}_{rbp}.pdf")))
