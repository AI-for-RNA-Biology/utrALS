# Notebooks and scripts

Which notebook makes which figures, what each notebook reads, and what the `utrals` package and the scripts do. How
the results read by the notebooks are generated: [`data_provenance.md`](data_provenance.md); how to regenerate them:
[`../pipeline/README.md`](../pipeline/README.md).

**Outputs.** Figures are written to `Figures/<notebook>/` (git-ignored; the `Figure_m6A_*` notebooks share
`Figures/Figure_m6A/`), tables to `Tables/<notebook>/` (versioned). `utrals.figures.save_figure` and `save_table` list
every saved file, with its legend or column description, in `Figures/README.md` and `Tables/README.md`.

**Data locations.** The notebooks import their shared code from the `utrals` package and read the data outside the
repository through `utrals.paths` (locations of `paths.env`).

**Reproducibility.** The notebooks fix their random seeds (`figures.fix_random_seeds`), save PDFs without a creation
date and contain nothing fetched live: what they read from the web is archived in `cache/` (STRING v12.5 responses,
Enrichr libraries, QuickGO annotations; delete a file to query again). With the package versions of `envs/`, a rerun
gives the same figures. `cache/reference_figures/` holds a GO bubble plot made from the live STRING website, whose GO
terms `Figure_3UTR_CTRL_vs_VCP` plots from the archived STRING v12.5 enrichment.

**Running the notebooks.** Kernel `utrals` for all, except `Figure_m6A_lncRNA_switch` (`bioinfo`). Run
`Figure_regulators_lengthening` before the `Figure_m6A_*` notebooks (§3). When running them headless, do not set
`MPLBACKEND=Agg` (the plots would not be stored in the notebooks), and set `R_HOME` for `Figure_3UTR_CTRL_vs_VCP` if R
is not found (see `envs/README.md`).

---

## 1. Figure notebooks

Named after their topic, not after a paper figure number. The "Contents" cell at the top of each notebook lists every
file it saves.

| Topic | Notebook(s) | Tables |
|---|---|---|
| 3′ UTR isoform usage over MN differentiation, nuclear vs cytoplasmic, CTRL vs VCP; RNA export | `Figure_3UTR_dynamics_over_time` | -- |
| 3′ UTR shifts in VCP vs CTRL MNs and in independent datasets (VCP-KI ± ML240, SFPQ / TDP-43 / FUS KD, ENCODE KD, astrocytes); pathways | `Figure_3UTR_CTRL_vs_VCP` | yes |
| Regulators selected at a single time point (Nuclear CTRL, D35) | `Figure_global_regulators_NuclearWT_D35` | yes |
| SVD of per-patient RBP activity over differentiation | `Figure_RBP_SVD` | -- |
| Regulators of lengthening (significant D22 → D35 change, Nuclear CTRL) | `Figure_regulators_lengthening` | yes |
| m6A, nucleotide content, m6A switches | `Figure_m6A_RBP_correlations`, `Figure_m6A_switch_regions_summary`, `Figure_m6A_lncRNA_switch`, `Figure_m6A_loading`, `Figure_m6A_nucleotide_content` | `Figure_m6A_nucleotide_content` |

---

## 2. What each notebook does and reads

Results paths are relative to the project folder (`PROJECT_DIR`); the pipeline step that writes each one is given in
brackets (MN_diff unless stated otherwise). External inputs are described in `data_provenance.md` §5.

| Notebook | What it does | Reads |
|---|---|---|
| `Figure_3UTR_dynamics_over_time` | nuclear vs cytoplasmic 3′ UTR isoform usage over MN differentiation, CTRL vs VCP; RNA export pathway genes over time | `MN_diff/results/3UTR_analysis/table_of_pairs.csv` [`bulkapa`]; `MN_diff/results/genes_preprocessing/df_genes_filtered_normalized.csv` [`genes_preprocessing`]; external: proteomics, RNA quantification, GO:0006406 genes (`cache/QuickGO/`) |
| `Figure_3UTR_CTRL_vs_VCP` | VCP-dependent loss of nuclear 3′ UTR lengthening; VCP-KI ± ML240, SFPQ / TDP-43 / FUS KD, ENCODE KD, astrocytes; STRING pathway enrichment | `table_of_pairs.csv` (and `shifts/`) of MN_diff, VCP-KI-Ml240, SFPQ_TDP_FUS_KD, Astrocytes and every `ENCODE-KD_*/<RBP>/` [`bulkapa`]; `VCP-KI-Ml240/results/genes_preprocessing/df_genes_filtered_normalized.csv` [`genes_preprocessing`] |
| `Figure_global_regulators_NuclearWT_D35` | top 15 regulators per region at a single time point (Nuclear CTRL, D35): STRING networks, pathway enrichment, expression, D22 → D35 change and time course, condensate propensity, VCP-KI, nuclear retention vs APA | `MN_diff/results/3UTR_analysis/CLIP/analysis/welch_{proximal,distal}.csv` [`clip_merge_welch`]; `MN_diff/results/RNA_modifs/<mod>/welch/analysis/welch_*.csv` [`rnamod_merge_welch`]; `genes_preprocessing`; `VCP-KI-Ml240/results/3UTR_analysis/CLIP/*.pkl` [`clip_integration`]; external: proteomics, AxonLoc |
| `Figure_RBP_SVD` | SVD of the per-patient CLIP Welch profiles per region and fraction | `MN_diff/results/3UTR_analysis/CLIP_per_patient/analysis/welch_*.csv` [`clip_per_patient_merge_welch`] |
| `Figure_regulators_lengthening` | RBP regulators of 3′ UTR lengthening (significant D22 → D35 change); VCP-KI and SFPQ / TDP-43 / FUS KD regulators; lncRNA cross-links; writes the candidate lists (§3) | as `Figure_global_regulators_NuclearWT_D35`, plus `MN_diff/results/paraspeckles/<GENE>/<GENE>_cross_links.csv` [`paraspeckles_*`], the m6A pickles of VCP-KI-Ml240 [`m6a_integration`], the CLIP and m6A pickles of SFPQ_TDP_FUS_KD [`clip_integration`, `m6a_integration`] |
| `Figure_m6A_RBP_correlations` | correlations between the CLIP and the m6A Welch profiles around the PAS | CLIP and RNA-modification `welch_*.csv` (as above); `MN_diff/results/RNA_modifs/<mod>/welch/corr_RBP_<mod>_welch_dist_min_150.csv` (§3); candidate lists (§3) |
| `Figure_m6A_switch_regions_summary` | m6A structural-switch competence and CLIP enrichment in the regions around the PAS | `MN_diff/results/RNA_modifs/m6A/regions/{comparison,clip_all_datasets}/` [`m6a_regions*`]; candidate lists |
| `Figure_m6A_lncRNA_switch` | m6A switches and CLIP preference along MALAT1, NEAT1_2, NORAD, XIST, TUG1 | `MN_diff/results/RNA_modifs/m6A/lncRNA_switch/` [`lncrna_*`]; candidate lists |
| `Figure_m6A_loading` | m6A density around the PAS and stop codon, normalised by A content; expression of the m6A machinery over time | `table_of_pairs.csv`, `shifts/` [`bulkapa`]; `genes_preprocessing`; external: RMBase m6A sites, 3′ UTR GTFs, hg38 genome |
| `Figure_m6A_nucleotide_content` | A/C/G/U and dinucleotide content around the PAS for the shift groups | `table_of_pairs.csv`, `shifts/` [`bulkapa`]; external: 3′ UTR GTF, hg38 genome |

---

## 3. Notebook outputs read by other analyses

| Notebook | Writes | Used by |
|---|---|---|
| `Figure_regulators_lengthening` | `Tables/Figure_regulators_lengthening/regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_<region>_APA_vs_nuclear_retention_class.csv` (`<region>` = `region1_upstream_pPAS`, `region2_downstream_pPAS`, `region3_upstream_dPAS`) | candidate RBPs of the m6A-switch pipeline (`m6a_regions` step) and the `Figure_m6A_*` notebooks |
| `Figure_m6A_RBP_correlations` | `MN_diff/results/RNA_modifs/<mod>/welch/corr_RBP_<mod>_welch_dist_min_150.{csv,xlsx}` (correlation of every RBP's Welch profile with the modification's; its cell is commented out, uncomment it to recompute) and `corr_RBP_m6a_welch_piled_*.csv` | the notebook itself |

---

## 4. The `utrals` package (`src/utrals/`)

| Module | Content | Used by |
|---|---|---|
| `paths` | data locations of `paths.env` and the repository folders (`REPO_ROOT`, `FIGURES_DIR`, `TABLES_DIR`, `CACHE_DIR`) | everything |
| `figures` | `save_figure`, `save_table` (and the `Figures/README.md`, `Tables/README.md` catalogues), `fix_random_seeds`, `lengthening_regulators_csv` | notebooks |
| `rbp_names` | CLIP dataset / RBP name parsing | notebooks, `welch` |
| `stats` | `find_elbow`, `get_outliers` | notebooks |
| `string_api`, `string_analysis` | STRING API requests (v12.5, archived in `cache/`), networks and enrichment bubble plots | notebooks |
| `enrichr` | Enrichr gene-set libraries (archived in `cache/Enrichr/`) | `Figure_3UTR_CTRL_vs_VCP` |
| `expression`, `go_expression` | gene and GO-term expression over the MN differentiation; WT vs VCP heatmaps | notebooks |
| `welch` | Welch profile figures: `panels` (3 × 3 cm panels), `fisher` (Fisher enrichment heatmaps), `common`, and one module per layout: `time_course` (MN_diff), `knockdown` (SFPQ_TDP_FUS_KD), `treatment_grid` (VCP-KI-Ml240) | `scripts/clip/clip_analysis.py` |
| `m6a_switch` | m6A structural-switch analysis (`competence`, `competence_gu_motif`, `group_stats`, `clip_enrichment_*`, `compare_regions`, `export_pairs_table`, `lncrna`; needs pysam and ViennaRNA) | `scripts/m6a_switch/`, `Figure_m6A_lncRNA_switch` |

---

## 5. Scripts (`scripts/`, run by the pipeline)

| Script | Purpose | Pipeline step |
|---|---|---|
| `counts/sum_UTR_counts_across_replicates.py`, `counts/sum_genes_counts_across_replicates.py` | sum the 3 sequencing runs of each MN_diff sample | `sum_utr_counts`, `sum_gene_counts` |
| `clip/clip_analysis.py --layout {time_course,knockdown,treatment_grid}` | Welch profiles of every CLIP dataset and RNA modification, and Fisher heatmaps, with the layout of the dataset | `clip_analysis`, `rnamod_analysis`, `m6a_analysis` |
| `clip/merge_welch.py`, `clip/merge_welch_per_patient.py` | Welch tables `analysis/welch_{proximal,distal}.csv` from the integration pickles | `clip_merge_welch`, `rnamod_merge_welch`, `clip_per_patient_merge_welch` |
| `clip/gene_site_density.py` | CLIP cross-links in 200-nt windows along a gene, and `<GENE>_cross_links.csv` | `paraspeckles_*` |
| `m6a_switch/run_m6a_switch_pipeline.sh` | m6A switch pipeline for one region (the scripts below, `bioinfo` environment) | `m6a_regions` |
| `m6a_switch/m6a_switch_competence_analysis.py`, `m6a_switch_competence_guu_motif_analysis.py` | RNAfold of the windows; switch-competent m6A sites (U-tract or GU motif) | via `run_m6a_switch_pipeline.sh` |
| `m6a_switch/switch_competence_group_stats.py` | nuclear-lengthened vs other statistics | via `run_m6a_switch_pipeline.sh` |
| `m6a_switch/clip_u_tract_switch_competent_enrichment.py`, `clip_guu_motif_switch_competent_enrichment.py` | CLIP enrichment of the candidate RBPs at switch-linked elements | via `run_m6a_switch_pipeline.sh` |
| `m6a_switch/clip_enrichment_all_datasets.py` | the same test for all CLIP datasets | `m6a_regions_clip_*` |
| `m6a_switch/compare_switch_regions.py`, `export_pairs_switch_table.py` | comparison of the regions; per-pair switch table | `m6a_regions_compare` |
| `m6a_switch/lncrna_m6a_switch_clip.py` | switch competence and CLIP preference along whole lncRNAs | `lncrna_*` |

Each `scripts/m6a_switch/*.py` script is the command-line entry point of the `utrals.m6a_switch` module of the same
topic.
