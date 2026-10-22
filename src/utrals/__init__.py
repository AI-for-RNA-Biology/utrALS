"""utrALS: code behind the paper figures on RBP binding and 3' UTR isoform usage in ALS motor neurons.

Modules
- paths: locations of the data outside the repository (paths.env) and of the repository folders
- figures: save_figure / save_table (Figures/README.md, Tables/README.md catalogues), fix_random_seeds
- rbp_names: CLIP dataset / RBP name parsing; stats: elbow and outlier helpers
- string_api: STRING API requests pinned to v12.5, archived in cache/; string_analysis: STRING networks and
  enrichment plots; enrichr: cached Enrichr gene-set libraries
- expression, go_expression: gene and GO-term expression over the MN differentiation time course
- welch: Welch profile figures of the integrated CLIP / RNA-modification results
- m6a_switch: m6A structural-switch analysis
"""
