import pandas as pd
import os
from utrals.paths import PROJECT_DIR as _PROJECT_DIR  # data locations, see paths.env
from cliptools.analysis import load_and_combine_integrated_CLIP

# Parse column names to extract day, compartment, and mutation
# Columns are like: -log10(pval)*dPUD_C_D3_Cytoplasmic_VCP or -log10(pval)*dPUD_M_D3_Cytoplasmic_WT
def parse_column_name(col):
    # Remove prefix '-log10(pval)*dPUD_'
    col_clean = col.replace('-log10(pval)*dPUD_', '')
    parts = col_clean.split('_')
    
    # parts[0] is C or M (compartment indicator)
    # parts[1] is Day (D0, D3, D7, etc)
    # parts[2] is Compartment (Cytoplasmic or Nuclear)
    # parts[3] is Mutation (VCP or WT)
    
    return {
        'Day': parts[0] if len(parts) > 0 else None,
        'Compartment': parts[1] if len(parts) > 1 else None,
        'VCP': parts[2] if len(parts) > 2 else None
    }
    
if __name__ == "__main__":
    import sys

    # Usage: python merge_welch.py [FOLDER_OF_INTEGRATED_PKL]   (writes <folder>/analysis/welch_{proximal,distal}.csv)
    # Used for MN_diff/results/3UTR_analysis/CLIP and MN_diff/results/RNA_modifs/m6A/welch (the default).
    if len(sys.argv) > 1:
        path_to_processed_RBPs = sys.argv[1]
    else:
        path_to_processed_RBPs = f"{_PROJECT_DIR}/MN_diff/results/RNA_modifs/m6A/welch"

    if not os.path.exists(os.path.join(path_to_processed_RBPs, "analysis")):
        os.makedirs(os.path.join(path_to_processed_RBPs, "analysis"))
        
    for isoform in ["proximal", "distal"]:
        print(f"Processing isoform: {isoform}")
        df_results = load_and_combine_integrated_CLIP(path_to_processed_RBPs, prefix=isoform)
        dfs_proxs = []
        for clip in df_results.keys():
            print(f"Processing CLIP: {clip}")
            d = df_results[clip][[col for col in df_results[clip].columns if col.startswith('-log10(pval)*dPUD_D')]]
            d.index = [-1000+i*50 for i in range(d.shape[0])]
            
            print(d)
            # Transform d from wide to long format
            d_long = d.reset_index().melt(id_vars='index', var_name='column_name', value_name='mean_score_lg10')
            d_long.columns = ['index', 'column_name', 'mean_score_lg10']



            # Apply parsing to each column name
            parsed_info = d_long['column_name'].apply(parse_column_name)
            d_long = d_long.assign(**{col: parsed_info.apply(lambda x: x[col]) for col in ['Day', 'Compartment', 'VCP']})
            print(d_long.head())
            # Clean up unnecessary columns
            d_long = d_long[['index', 'Day', 'Compartment', 'VCP', 'mean_score_lg10']]
            d_long = d_long.rename(columns={'index': 'dist'})
            d_long["Isoform"] = f"{isoform.capitalize()} 3'UTR Isoform"
            d_long["CLIP"] = clip

            dfs_proxs.append(d_long)
    
        df_tot_prox = pd.concat(dfs_proxs)
        df_tot_prox.to_csv(os.path.join(path_to_processed_RBPs, "analysis", f"welch_{isoform}.csv"))