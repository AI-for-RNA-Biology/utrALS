"""Welch profiles on a fraction x (stop codon / proximal / distal) grid, one line per genotype x drug treatment
(VCP-KI-Ml240).

Layout `treatment_grid` of scripts/clip/clip_analysis.py.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from .panels import CM, RC, LINEWIDTH, panel_grid, figure_title
from .common import setup_output_directories, load_integrated, save, for_each


def split_treatment(treatment):
    """'Ml240_VCP_KI_Cyto' -> ('Ml240', 'VCP_KI')."""
    parts = treatment.split('_')
    return parts[0], '_'.join(parts[1:-1])


def clipped_profiles(df_results, treatments, sign=1):
    """
    -log10(pval)*dPUD per treatment. pvalue == 0 gives +-inf: cap it at the most extreme
    finite value observed so a degenerate bin doesn't blow up the y-axis.
    """
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in treatments]
    raw = df_results[cols].apply(pd.to_numeric, errors='coerce')
    finite = raw.replace([np.inf, -np.inf], np.nan)
    return raw.replace([np.inf, -np.inf], [finite.max().max(), finite.min().min()]) * sign


def plot_welch_profiles(df_proximal, df_distal, RBP, config, df_stop_codon=None):
    """
    Welch profiles (lines only, panels.PANEL_CM panels): rows = compartments,
    columns = stop codon (if `df_stop_codon` is given) / proximal / distal.
    Save with bbox_inches='tight'.
    """
    treatments = config['treatments']
    compartments = config['compartments']
    palette = config['palette_conditions']
    nt_before = config.get('nt_before', -1000)
    length_region = config.get('length_region', 50)

    # Distal profiles are plotted as -value (same convention as the MN_diff Welch profiles);
    # stop codon profiles keep their sign, like the proximal ones (positive = associated with proximal usage)
    columns = [(df_proximal, 'pPAS', 1, 'Proximal', nt_before, length_region),
               (df_distal, 'dPAS', -1, 'Distal', nt_before, length_region)]
    if df_stop_codon is not None:
        columns.insert(0, (df_stop_codon, 'stop codon', 1, 'Stop codon',
                           config.get('stop_codon_nt_before', -100), config.get('stop_codon_length_region', 25)))

    with plt.rc_context(RC):
        fig, axes = panel_grid(len(compartments), len(columns), right_cm=4.0)
        for col, (df, anchor, sign, name, col_nt_before, col_length_region) in enumerate(columns):
            clipped = clipped_profiles(df, treatments, sign)
            ylim = (clipped.min().min() - 2, clipped.max().max() + 2)
            distances = [col_nt_before + i * col_length_region for i in range(df.shape[0])]
            for row, (compartment, key) in enumerate(compartments.items()):
                ax = axes[row, col]
                for treatment in treatments:
                    if not treatment.endswith(key):
                        continue
                    drug, genotype = split_treatment(treatment)
                    ax.plot(distances, clipped[f'-log10(pval)*dPUD_{treatment}'],
                            label=f'{genotype} {drug}', color=palette[f'{drug}_{genotype}'], linewidth=LINEWIDTH)
                ax.axhline(0, color='grey', linestyle=':', linewidth=0.5)
                ax.axvline(0, color='grey', linestyle=':', linewidth=0.5)
                if row == len(compartments) - 1:
                    ax.set_xlabel(f'Distance to {anchor} (nt)')
                if col == 0:
                    ax.set_ylabel('-log10(pval)*dPUD')
                ax.set_title(f'{compartment} - {name}')
                ax.set_ylim(ylim)
                sns.despine(ax=ax)

        handles, labels = axes[0, 0].get_legend_handles_labels()
        right, middle = axes[0, -1].get_position().x1, (axes[0, -1].get_position().y1 + axes[-1, -1].get_position().y0) / 2
        fig.legend(handles, labels, loc='center left', frameon=False,
                   bbox_to_anchor=(right + 0.4 * CM / fig.get_figwidth(), middle))
        figure_title(fig, axes, RBP, offset_cm=0.5)
    return fig



def run(config, welch_only=False):
    """Welch profiles of one analysis config (see pipeline/datasets/VCP-KI-Ml240/configs/*analysis*.yaml); only the
    columns needed for the plots are loaded."""
    out = setup_output_directories(config['output_directory'], {'welch': 'welch_profiles'})
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in config['treatments']]
    proximals, distals, stop_codon = load_integrated(config, cols=cols, RBP_list=config.get('RBP_list'))
    rbps_common = sorted(rbp for rbp in proximals if rbp in distals)
    missing = set(proximals) ^ set(distals)
    if missing:
        print(f"  Skipping {len(missing)} RBPs without both proximal and distal results: {sorted(missing)}")
    for_each(rbps_common, "Welch profiles (stop codon + proximal + distal)", lambda rbp: save(
        out['welch'] / f"welch_profile_{rbp}.pdf",
        plot_welch_profiles(proximals[rbp], distals[rbp], rbp, config, df_stop_codon=stop_codon.get(rbp))))
