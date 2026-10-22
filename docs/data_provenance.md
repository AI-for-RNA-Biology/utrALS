# Data provenance: how the results read by the figures are generated

The project folder (`PROJECT_DIR` in `paths.env`) holds the processed data and results of the ALS 3′ UTR / RBP
compartmentalisation project: alternative polyadenylation (APA), CLIP-based RBP binding, RNA modifications and gene
expression, measured in nuclear and cytoplasmic fractions across several ALS-relevant datasets. Paths below are relative
to that folder unless stated otherwise; `$ANNOT_DIR`, `$CLIP_DIR`, ... are the locations defined in `paths.env`.

This document covers the results read by the paper-figure notebooks (`notebooks/`); which notebook reads which result
is listed in [`notebooks_and_scripts.md`](notebooks_and_scripts.md) §2. Every dataset is regenerated with
`pipeline/datasets/<dataset>/run_all.sh` (see [`../pipeline/README.md`](../pipeline/README.md)).

**Project owner**: Lisa Fournier, lisa.fournier@unibe.ch

---

## 1. Datasets

| Folder | Dataset | Raw data source |
|---|---|---|
| `MN_diff/` | iPSC → motor neuron differentiation, VCP-mutant vs control, nuclear/cytoplasmic fractions, D0–D35 | GEO GSE152983 / SRA PRJNA641118 (Tyzack et al. 2021, PMID 33693641) |
| `Astrocytes/` | iPSC-derived astrocytes, VCP vs control, nuclear/cytoplasmic | SRA PRJNA672040 |
| `VCP-KI-Ml240/` | VCP knock-in / VCP / TDP-43 / WT lines ± ML240 (VCP inhibitor), nuclear/cytoplasmic | SRA PRJNA883641 |
| `SFPQ_TDP_FUS_KD/` | siRNA knockdown of SFPQ, TDP-43 and FUS (+ untreated, scramble), nuclear/cytoplasmic | In-house, soon released on GEO |
| `ENCODE-KD_K562/`, `ENCODE-KD_HepG2/` | ENCODE shRNA knockdowns of single RBPs vs control | ENCODE (via SRA) |
| `AxonLoc/` | `Supplementary_tables.csv` (Excel workbook): condensate propensity z-scores of RBPs, read by two notebooks (§5) | AxonLoc supplementary tables |

References:
- **Astrocytes** (PRJNA672040): Ziff, O. J. et al. Meta-analysis of human and mouse ALS astrocytes reveals multi-omic
  signatures of inflammatory reactive states. Genome Res. 32, 71–84 (2022).
- **VCP-KI-Ml240** (PRJNA883641): Ziff, O. J. et al. Nucleocytoplasmic mRNA redistribution accompanies RNA binding
  protein mislocalization in ALS motor neurons and is restored by VCP ATPase inhibition. Neuron 111, 3011-3027.e7 (2023).
- **SFPQ_TDP_FUS_KD**: will be released on GEO soon.
- **AxonLoc supplementary tables**: Luisier, R., Andreassi, C., Fournier, L. & Riccio, A. The predicted RNA-binding
  protein regulome of axonal mRNAs. Genome Res. 33, 1497–1512 (2023).

Every dataset folder (except ENCODE) is split into `data/` (alignment and quantification) and `results/`
(`genes_preprocessing/`, `3UTR_analysis/` for APA, `RNA_modifs/` for RNA modifications, and dataset-specific folders).

Software: the conda environments of `envs/` and the AI-for-RNA-Biology packages at release `v0.2.0` (see
[`../envs/README.md`](../envs/README.md)).

---

## 2. Common pipeline

All bulk RNA-seq datasets go through the same steps; §3 lists the differences. The pipeline step that runs each one is
given in brackets.

**Step 1: Alignment (STAR)** → `data/STAR_hg38/` [`star`]
```
STAR --genomeDir genome_index_99 --readFilesIn R1 R2 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
```

**Step 2: Transcript quantification (kallisto 0.46.1)** → `data/hg38_kallisto/` [`kallisto`]
```
kallisto quant --index hg38_kallisto_v13.idx --rf-stranded -b 100 R1 R2
```

**Step 3: 3′ UTR isoform counting (featureCounts)** → `data/myUTRs_hg38/` [`featurecounts`, `aggregate`]
1. RSeQC `infer_experiment.py` on the BAM sets the strandedness: stranded if the dominant fraction is ≥ 0.55 and at
   least 2× the other fraction, or if it is ≥ 0.60; unstranded otherwise.
2. `featureCounts -a hg38_NoOverlap.gtf -g transcript_id -f -O -M -s <strandedness> -p`
3. All samples are merged into `df_counts.csv` (`bulkanalysis/scripts/aggregate_featureCounts_output.py`).

**Step 4: Gene-level preprocessing** → `results/genes_preprocessing/` [`genes_preprocessing`]
`bulkanalysis/scripts/genes_preprocessing.py`, configured by the `genes_preprocessing.yaml` stored in the output folder:
1. merge the kallisto outputs and sum the transcripts per protein-coding gene (GENCODE v44 names, NCBI gene_info);
2. per sample, fit a Gaussian mixture on log2(counts + 1) to set an expression threshold;
3. keep the genes above the threshold in ≥ `pct_in_treatment` of the samples of at least one condition;
4. quantile-normalise.

Outputs: `df_genes_filtered_normalized.csv`, `DGE/` (inputs of a differential expression analysis), `Figures/`,
`message.log` (stdout of the run).

**Step 5: 3′ UTR / APA analysis** → `results/3UTR_analysis/` [`bulkapa`]
`apatools/scripts/bulkapa.py`, configured by the `bulkapa.yaml` stored in the output folder:
1. correct the counts of short proximal isoforms (< 300 nt) and remove isoforms with outlier lengths;
2. keep the multi-UTR isoforms of the genes kept in Step 4, and flag genes with abnormal proximal abundance
   (`insane_genes_proximal_abundance.csv`);
3. apply the bulk 3′-bias correction;
4. filter isoforms per sample (Gaussian mixture threshold, `pct_in_treatment`) and quantile-normalise
   (`df_counts_filtered_corrected_normalized.csv`);
5. keep the shortest and longest isoform of each gene as its proximal / distal pair (`table_of_pairs_simplified.csv`);
6. for each pair of conditions, compute the PUD (proximal usage), RUD and DUD scores and run Fisher tests. A shift is
   significant when FDR < 0.01, |ΔPUD| > 0.15, and the RUD (1) and p-value (0.05) criteria are met → `table_of_pairs.csv`.

Outputs: `table_of_pairs.csv`; `shifts/pairs_{proximals,distals}_<A>_vs_<B>.csv` (pairs shifting to proximal usage,
i.e. shortening, or distal usage, i.e. lengthening, in B) and the matching `genes_*.csv`; `genes_utrs/`; `Figures/`;
`message.log` (including the number of significant pairs per comparison).

**Step 6: CLIP × APA integration** → `results/3UTR_analysis/CLIP/*.pkl` [`clip_integration`]
`apatools/scripts/clip_integration.py`, one run per CLIP file and isoform type (configs `CLIP/config_{eCLIP,iCLIP}.yaml`,
which differ only in how the RBP name is parsed from the file name). 335 CLIP datasets: 140 eCLIP HepG2, 177 eCLIP
K562, 18 iCLIP (selection in `list_clip_files`, `pipeline/lib.sh`). For each RBP and isoform type:
1. tile the region around the PAS into windows and flag each isoform as cross-linked or not in each window;
2. Welch tests on the PUD of cross-linked vs non-cross-linked pairs, per condition and between conditions;
3. Fisher tests of the association between cross-linking and the significant shift lists.

Output: `<RBP>_<rep>_<cellline>_results_{proximal,distal}.pkl`.

The same integration is run around the stop codon (`--iso_type stop_codon`, configs `CLIP/config_stop_codon_*.yaml`;
step `clip_stop_codon_integration`): windows from −100 to +100 nt (25-nt regions) around the 3′ UTR start (first
nucleotide after the stop codon). The proximal and distal isoforms of a pair share that start, so each pair is counted
once, under its distal isoform → `*_results_stop_codon.pkl`.

**Step 7: Welch profiles** → `results/3UTR_analysis/CLIP/analysis*/`, `results/RNA_modifs/<mod>/welch/analysis*/`
[`clip_analysis`, `rnamod_analysis`, `m6a_analysis`]
`scripts/clip/clip_analysis.py --layout <layout>` (code in `utrals.welch`; layout `time_course` for MN_diff,
`treatment_grid` for VCP-KI-Ml240, `knockdown` for SFPQ_TDP_FUS_KD). Welch profiles from −1000 to +200 nt around the
PAS in 50-nt regions, plus the stop codon profile (sign as the proximal one: positive = associated with proximal usage).
These are the Welch profiles shown in the paper. With `--welch-only` (MN_diff CLIP) only the profiles over time are
drawn; without it, the by-region profiles and the Fisher outputs are drawn too.

**Step 8: Welch tables (MN_diff)** → `analysis/welch_{proximal,distal}.csv` [`clip_merge_welch`, `rnamod_merge_welch`]
`scripts/clip/merge_welch.py` gathers the Welch scores (−log10(p) × ΔPUD) of every CLIP dataset (or RNA modification),
condition and 50-nt window into one table per isoform type. These tables are read by the notebooks.

**Step 9: RNA modifications × APA** → `results/RNA_modifs/<mod>/welch/` [`rnamod_integration`, `m6a_integration`]
Step 6 with the RMBase site file of each mark as the "CLIP file" (e.g. `RMBase/human.hg38.m6A.gff`), at the PAS and the
stop codon, followed by Steps 7 and 8.

---

## 3. Datasets

### 3.1 `MN_diff/` (main dataset)

`pipeline/datasets/MN_diff/run_all.sh`; samples in `samples.tsv`.

**Design**
- 288 SRA runs (SRR12068258–SRR12068545) → 96 samples: 8 lines × 6 days (D0, D3, D7, D14, D22, D35) × 2 fractions
  (Nuclear, Cytoplasmic).
- Lines: CTRL1, CTRL3, CTRL4, CTRL5 (wild type); MUT1, MUT2 (VCP R155C); MUT3, MUT4 (VCP R191Q).
- The 3 sequencing runs of each sample are summed, not treated as replicates
  (`scripts/counts/sum_{UTR,genes}_counts_across_replicates.py`; steps `sum_utr_counts`, `sum_gene_counts`).
- Conditions group the 4 lines of a genotype: `D<day>_<Fraction>_<Wildtype|VCP>` (24 conditions).

**`results/`**

| Folder | How it is generated |
|---|---|
| `genes_preprocessing/` | Step 4, 24 conditions |
| `3UTR_analysis/` | Step 5, 60 condition pairs (`treatment_pairs` in `bulkapa.yaml`): consecutive and D0-vs-later days within each fraction and genotype, Nuclear vs Cytoplasmic, Wildtype vs VCP at each day |
| `3UTR_analysis/CLIP/` | Steps 6–8, at the PAS and the stop codon |
| `3UTR_analysis/CLIP_per_patient/` | Steps 6 and 8 with the PUD computed per line instead of per condition (92 of the 96 samples, Welch tests only; `scripts/clip/merge_welch_per_patient.py`; steps `clip_per_patient_integration`, `clip_per_patient_merge_welch`). Used by `Figure_RBP_SVD` |
| `paraspeckles/<GENE>/` (MALAT1, NEAT1, NORAD, SFPQ, XIST) | `scripts/clip/gene_site_density.py` (steps `paraspeckles_*`): CLIP cross-links in 200-nt windows along the gene (GENCODE v44) and the per-CLIP totals `<GENE>_cross_links.csv` |
| `RNA_modifs/{m6A,m1A,m5C,m7G,Nm,Pseudo,RNA-editing}/welch/` | Step 9 with the RMBase sites of each mark |
| `RNA_modifs/m6A/regions/<region>/` | m6A structural-switch pipeline (below), step `m6a_regions` |
| `RNA_modifs/m6A/regions/clip_all_datasets/` | the CLIP enrichment test of the switch pipeline for every CLIP dataset, to rank the candidates against all RBPs (steps `m6a_regions_clip_*`) |
| `RNA_modifs/m6A/regions/comparison/` | comparison of the regions and per-pair switch table (step `m6a_regions_compare`) |
| `RNA_modifs/m6A/lncRNA_switch/` | m6A switch competence and CLIP preference along the mature MALAT1, NEAT1_2, NORAD, XIST and TUG1 transcripts (`scripts/m6a_switch/lncrna_m6a_switch_clip.py`, steps `lncrna_*`) |

**m6A structural-switch pipeline** (`scripts/m6a_switch/run_m6a_switch_pipeline.sh`), run for 4 regions: 100 nt on
each side of the 3′ UTR start (no CLIP test), 350 nt upstream of the proximal PAS, 200 nt downstream of the proximal
PAS, 350 nt upstream of the distal PAS. Per region:
1. fold the window with RNAfold (ViennaRNA); an RMBase m6A site is switch-competent if it pairs (p ≥ 0.5) with a U in a
   run of ≥ 3 U's → `switch_competence_u_tract/`;
2. the same with UGU/UUG/GUG/GUU motifs → `switch_competence_gu_motif/`;
3. statistics of nuclear-lengthened vs other pairs → `group_stats.csv`;
4. CLIP enrichment of the candidate RBPs at switch-linked U-tracts and GU motifs → `clip_u_tract_enrichment/`,
   `clip_gu_motif_enrichment/`.

The candidate RBPs are the regulators of lengthening written by `notebooks/Figure_regulators_lengthening.ipynb`
(`Tables/Figure_regulators_lengthening/regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_<region>_APA_vs_nuclear_retention_class.csv`),
so that notebook runs before this step. Each step writes a `params.json` with its inputs and package versions.

**Proteomics** of the same differentiation (GEO GSE152983, column-median normalised): `$PROTEOMICS_CSV`, read by the
notebooks; not produced by the pipeline.

### 3.2 `Astrocytes/`

`pipeline/datasets/Astrocytes/run_all.sh`. 8 runs (SRR12899403–SRR12899410) in 4 conditions (`CTRL_CYTO`, `VCP_CYTO`,
`CTRL_NUC`, `VCP_NUC`). Steps 1–5, all pairwise comparisons.

### 3.3 `VCP-KI-Ml240/`

`pipeline/datasets/VCP-KI-Ml240/run_all.sh`. 47 runs (SRR21683985–SRR21684031) in 16 conditions:
`{Ml240, Untreated}` × `{TDP43, VCP_KI, VCP, WT}` × `{Cyto, Nuc}`. Steps 1–5 (selected `treatment_pairs`), Steps 6–7
for CLIP and Step 9 for m6A, at the PAS and the stop codon.

### 3.4 `SFPQ_TDP_FUS_KD/`

`pipeline/datasets/SFPQ_TDP_FUS_KD/run_all.sh`. 30 samples, 3 replicates per condition, 10 conditions:
`{Untreated, Scramble, SFPQ_KD, TDP43_KD, FUS_KD}` × `{Cyto, Nuc}`; paired-end over 2 lanes. Steps 1–5, Steps 6–7 for
CLIP and Step 9 for m6A, at the PAS and the stop codon.

### 3.5 `ENCODE-KD_K562/` and `ENCODE-KD_HepG2/`

`pipeline/datasets/ENCODE-KD_{K562,HepG2}/run_all.sh` (`ENCODE_RBPS="TARDBP FUS ..."` runs only some RBPs). One folder
per knocked-down RBP, comparing `control` with `<RBP>-KD`; each cell line has one shared control. Analysed RBPs: 14 in
HepG2; FUS, SFPQ and TARDBP in K562. Strandedness is inferred per run with RSeQC for both kallisto and featureCounts; no
kallisto bootstraps. Steps 1–3 are written to `$EXTERNAL_DATA_DIR/{K562,HepG2}_shRNA/`; `<RBP>/genes_preprocessing/`
and `<RBP>/3UTR_analysis/` are Steps 4–5 (`pct_in_treatment` 0.5, `quantile_threshold` 0.05).

---

## 4. Reference data

Read by the pipeline: genome indices and annotations (`$ANNOT_DIR`), CLIP cross-link files (`$CLIP_DIR`: eCLIP HepG2
and K562, iCLIP) and RMBase RNA-modification sites (`$RMBASE_DIR`).

---

## 5. External inputs of the notebooks

Tables read directly by the notebooks, not produced by the pipeline (their locations are variables of `paths.env`):
- proteomics GSE152983 (Tyzack et al. 2021): `$PROTEOMICS_CSV` (`Figure_3UTR_dynamics_over_time`,
  `Figure_global_regulators_NuclearWT_D35`, `Figure_regulators_lengthening`);
- RNA quantification of MN_diff: `$DATABASE_DIR/PRJNA641118/RNA_quantification.csv` (`Figure_3UTR_dynamics_over_time`);
- AxonLoc supplementary tables, condensate propensity z-scores (6th sheet; an Excel workbook despite its `.csv` name):
  `$AXONLOC_CSV` (`Figure_global_regulators_NuclearWT_D35`, `Figure_regulators_lengthening`);
- RMBase m6A sites, 3′ UTR GTFs and the hg38 genome (`$ANNOT_DIR`) (`Figure_m6A_loading`, `Figure_m6A_nucleotide_content`).
