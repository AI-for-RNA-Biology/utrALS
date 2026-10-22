"""m6A structural-switch analysis (MN_diff): RNAfold of windows around the PAS, switch-competent m6A sites
(paired with a U-tract or a UG motif), CLIP enrichment of candidate RBPs at the switch-linked elements, comparison of
the 3' UTR regions, and the same along whole lncRNAs. Command-line entry points: scripts/m6a_switch/.

Modules (command-line script in brackets)
- competence (m6a_switch_competence_analysis.py), competence_gu_motif (m6a_switch_competence_guu_motif_analysis.py)
- group_stats (switch_competence_group_stats.py)
- clip_enrichment_u_tract (clip_u_tract_switch_competent_enrichment.py),
  clip_enrichment_gu_motif (clip_guu_motif_switch_competent_enrichment.py),
  clip_enrichment_all_datasets (clip_enrichment_all_datasets.py)
- compare_regions (compare_switch_regions.py), export_pairs_table (export_pairs_switch_table.py)
- lncrna (lncrna_m6a_switch_clip.py)
"""
