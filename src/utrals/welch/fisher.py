"""Fisher tests of the CLIP / APA integration: signed enrichment of cross-links in the shifted pairs, heatmaps."""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ..rbp_names import get_clip_name


def compute_signed_enrichment(fisher_results_series):
    """
    Compute signed enrichment score from fisher test results.
    
    Parameters:
    -----------
    fisher_results_series : pandas.Series
        Series where each element is a tuple (odds_ratio, pvalue)
        
    Returns:
    --------
    pandas.Series
        Series with signed enrichment scores: sign(odds_ratio - 1) * -log10(pvalue)
        Positive values: enrichment of crosslinks in shifted pairs
        Negative values: depletion of crosslinks in shifted pairs
    """
    signed_enrichment = []
    for result in fisher_results_series:
        odds_ratio, pvalue = result
        # Determine sign: +1 if enriched (odds_ratio > 1), -1 if depleted (odds_ratio < 1)
        sign = 1 if odds_ratio > 0 else -1
        # Compute signed enrichment score
        score = sign * -np.log10(pvalue)
        signed_enrichment.append(score)
    
    return pd.Series(signed_enrichment, index=fisher_results_series.index)



def time_course_title(pairs_name):
    """Heatmap title of an MN_diff pair, e.g. pairs_distals_D22_Cytoplasmic_VCP_vs_D35_Cytoplasmic_VCP."""
    parts = pairs_name.split('_')
    day1, compartment, mutation, day2 = parts[2], parts[3], parts[4], parts[6]
    return f"Fisher ({day2} vs {day1}) \n {compartment} {mutation}"


def treatment_title(pairs_name):
    """Heatmap title of a pair of treatments, e.g. pairs_distals_Scramble_Nuc_vs_TDP43_KD_Nuc."""
    full_name = pairs_name.split('als_')[1]
    treatment1 = full_name.split('_vs_')[0]
    treatment2 = full_name.split('_vs_')[1]
    return f"Fisher ({treatment2} vs {treatment1})"


def plot_fisher_heatmap(pairs_name,
                        df_results,
                        title,
                        top_n=50,
                        nb_nucleotides_upstream=-1000,
                        length_region=50,
                        save_path=None,
                        limits=None,
                        return_matrix=False,
                        merge_ptb=True):
    """
    Generate a heatmap of Fisher enrichment test results for RBPs.

    Parameters:
    -----------
    pairs_name : str
        Name of the pairs file (e.g., "pairs_distals_D22_Cytoplasmic_VCP_vs_D35_Cytoplasmic_VCP")
    df_results : dict
        Dictionary of DataFrames containing fisher test results (e.g., df_results_proximals or df_results_distals)
    title : str
        Plot title (time_course_title / treatment_title below)
    top_n : int, optional
        Number of top RBPs to display (default: 50)
    merge_ptb : bool
        Passed to rbp_names.get_clip_name: merge every PTB / PTBP1 CLIP into one "PTB" row
    save_path : str, optional
        Path to save the figure
    """
    distances = [nb_nucleotides_upstream+i*length_region for i in range(df_results[list(df_results.keys())[0]].shape[0])]

    reference_clip = next(iter(df_results))
    reference_index = compute_signed_enrichment(
        df_results[reference_clip][f"fisher_results_{pairs_name}"]
    ).index

    matrix = np.vstack([
        compute_signed_enrichment(df_results[clip][f"fisher_results_{pairs_name}"])
        .reindex(reference_index)
        .to_numpy()
        for clip in df_results.keys()
    ])

    df_matrix = pd.DataFrame(matrix, index=df_results.keys(), columns=distances)

    if limits is not None:
        cols_to_keep = [limits[0]+i*length_region for i in range((limits[1]-limits[0])//length_region)]
        df_matrix = df_matrix[cols_to_keep]

    # Sort by maximum value across rows (highest at top)
    df_matrix = df_matrix.loc[df_matrix.max(axis=1).sort_values(ascending=False).index]

    # Map RBP names
    df_matrix["RBP_name"] = df_matrix.index.map(lambda clip: get_clip_name(clip, merge_ptb=merge_ptb))
    # df_matrix["RBP_name"].replace("TARDBP", "TDP43", inplace=True)
    df_matrix_merged = df_matrix.groupby("RBP_name").mean()
    df_matrix_merged = df_matrix_merged.loc[df_matrix_merged.max(axis=1).sort_values(ascending=False).index]

    # Create mask for non-significant values (absolute value < -log10(0.05) ≈ 1.301)
    significance_threshold = -np.log10(0.05)
    mask = np.abs(df_matrix_merged.head(top_n)) < significance_threshold

    # Get the maximum absolute value across the entire dataframe
    max_abs_val = df_matrix_merged.head(top_n).abs().max().max()

    plt.figure(figsize=(8, 8))
    sns.heatmap(df_matrix_merged.head(top_n), cmap="RdBu_r", mask=mask, yticklabels=True, center=0,
                vmin=-max_abs_val, vmax=max_abs_val)
    plt.title(title)

    if save_path:
        plt.savefig(save_path, bbox_inches="tight", metadata={"CreationDate": None} if str(save_path).endswith(".pdf") else None)

    if return_matrix:
        return df_matrix_merged
