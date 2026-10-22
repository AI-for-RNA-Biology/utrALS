"""
Fixed-size panel grids shared by the Welch profile figures of all datasets
(scripts/clip/clip_analysis.py): every panel is PANEL_CM x PANEL_CM,
whatever the number of panels, so that the profiles of different datasets have the same scale.

Create the figure and its artists inside `plt.rc_context(RC)`, and save with
bbox_inches='tight': it crops the page to its content without rescaling the panels.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

CM = 1 / 2.54
PANEL_CM = 3.0
LINEWIDTH = 1.0
RC = {'font.size': 6, 'axes.titlesize': 6, 'axes.labelsize': 6, 'xtick.labelsize': 5,
      'ytick.labelsize': 5, 'legend.fontsize': 6, 'legend.title_fontsize': 6, 'figure.titlesize': 8,
      'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
      'xtick.major.size': 2.5, 'ytick.major.size': 2.5}


def panel_grid(nrows, ncols, panel_cm=PANEL_CM, wspace_cm=1.1, hspace_cm=1.0, block_size=None,
               block_gap_cm=0.8, left_cm=1.2, right_cm=0.5, top_cm=1.0, bottom_cm=1.0):
    """
    Figure with a nrows x ncols array of axes, each exactly `panel_cm` x `panel_cm`.
    Columns can be grouped in blocks of `block_size` columns separated by an extra `block_gap_cm`
    (see block_title). The margins and spaces are in cm and leave room for the labels and titles.
    """
    n_block_gaps = (ncols - 1) // block_size if block_size else 0
    width = left_cm + ncols * panel_cm + (ncols - 1) * wspace_cm + n_block_gaps * block_gap_cm + right_cm
    height = top_cm + nrows * panel_cm + (nrows - 1) * hspace_cm + bottom_cm
    fig = plt.figure(figsize=(width * CM, height * CM))
    axes = np.empty((nrows, ncols), dtype=object)
    for row in range(nrows):
        for col in range(ncols):
            x = left_cm + col * (panel_cm + wspace_cm) + (col // block_size if block_size else 0) * block_gap_cm
            y = height - top_cm - (row + 1) * panel_cm - row * hspace_cm
            axes[row, col] = fig.add_axes([x / width, y / height, panel_cm / width, panel_cm / height])
    return fig, axes


def block_title(fig, first_ax, last_ax, title, offset_cm=0.55, **kwargs):
    """Title centred above the panels from `first_ax` to `last_ax` (same row), above their own titles."""
    left, right = first_ax.get_position().x0, last_ax.get_position().x1
    top = first_ax.get_position().y1 + offset_cm * CM / fig.get_figheight()
    fig.text((left + right) / 2, top, title, ha='center', va='bottom', fontweight='bold', **kwargs)


def figure_title(fig, axes, title, offset_cm=1.0):
    """Figure title centred above the whole grid."""
    left, right = axes[0, 0].get_position().x0, axes[0, -1].get_position().x1
    top = axes[0, 0].get_position().y1 + offset_cm * CM / fig.get_figheight()
    fig.text((left + right) / 2, top, title, ha='center', va='bottom', fontweight='bold',
             fontsize=plt.rcParams['figure.titlesize'])


def single_welch_profile(df_results, RBP, treatments, palette, anchor='pPAS', sign=1,
                         nt_before=-1000, length_region=50):
    """
    One PANEL_CM x PANEL_CM Welch profile panel (lines only), one line per treatment,
    legend on the right. `sign` flips the profile (distal profiles are plotted as -value).
    Save with bbox_inches='tight'.
    """
    cols = [f'-log10(pval)*dPUD_{treatment}' for treatment in treatments]
    # pvalue == 0 gives +-inf: cap it at the most extreme finite value observed
    raw = df_results[RBP][cols].apply(pd.to_numeric, errors='coerce')
    finite = raw.replace([np.inf, -np.inf], np.nan)
    clipped = raw.replace([np.inf, -np.inf], [finite.max().max(), finite.min().min()]) * sign
    distances = [nt_before + i * length_region for i in range(df_results[RBP].shape[0])]

    with plt.rc_context(RC):
        fig, axes = panel_grid(1, 1, right_cm=4.0, top_cm=0.8)
        ax = axes[0, 0]
        for treatment in treatments:
            ax.plot(distances, clipped[f'-log10(pval)*dPUD_{treatment}'], label=treatment,
                    color=palette[treatment], linewidth=LINEWIDTH)
        ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
        ax.set_xlabel(f'Distance to {anchor} (nt)')
        ax.set_ylabel('-log10(pval)*dPUD')
        ax.set_title(RBP)
        ax.set_ylim(clipped.min().min() - 2, clipped.max().max() + 2)
        sns.despine(ax=ax)
        ax.legend(title='Treatments', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    return fig
