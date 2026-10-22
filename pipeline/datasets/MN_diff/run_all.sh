#!/bin/bash
# Regenerate PROJECT_DIR/MN_diff (data/ and results/) from the raw FASTQ files.
#
# Usage: bash run_all.sh [--from STEP] [--to STEP] [--only STEP,...] [--jobs N] [--dry-run]
#        bash run_all.sh --list-steps
# On SLURM: bash pipeline/slurm/submit.sh MN_diff [same step options]
#
# Only the results read by the paper-figure notebooks (notebooks/).
# Not included (see docs/data_provenance.md): the analyses done in notebooks (RNA_modifs/*/welch/corr_RBP_*).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/params.sh"
source "${PIPELINE_DIR}/runner.sh"

STEPS=(download star kallisto featurecounts aggregate sum_utr_counts sum_gene_counts
       genes_preprocessing bulkapa
       clip_integration clip_stop_codon_integration clip_analysis clip_merge_welch
       clip_per_patient_integration clip_per_patient_merge_welch
       paraspeckles_prepare paraspeckles_run paraspeckles_merge
       rnamod_integration rnamod_stop_codon_integration rnamod_analysis rnamod_merge_welch
       m6a_regions m6a_regions_clip_list m6a_regions_clip_run m6a_regions_clip_merge m6a_regions_compare
       lncrna_prepare lncrna_features lncrna_run lncrna_merge)

# ── samples ──────────────────────────────────────────────────────────────────
# runs: every sequencing run (3 per sample); the counts of a sample's runs are summed (sum_* steps)
runs() { awk -F'\t' 'NR > 1 && NF {n = split($2, r, ","); for (i = 1; i <= n; i++) print r[i]}' "${DATASET_DIR}/samples.tsv"; }
# groups COLUMN: print "<group> <run> <run> ..." for each value of a samples.tsv column (1 = sample, 7 = category)
groups() {
    awk -F'\t' -v c="$1" 'NR > 1 && NF { if (!($c in g)) order[++n] = $c; gsub(",", " ", $2); g[$c] = g[$c] " " $2 }
        END { for (i = 1; i <= n; i++) print order[i] g[order[i]] }' "${DATASET_DIR}/samples.tsv"
}
fq() { echo "${FASTQ_DIR}/$1_$2.fastq.gz"; }
bam() { echo "${DATA_DIR}/STAR_hg38/$1Aligned.sortedByCoord.out.bam"; }
py() { task run_script ENV_PY python "$@"; }          # project script in the utrals environment
py_bio() { task run_script ENV_BIO python "$@"; }     # project script in the bioinfo environment (m6A switch)

# ── raw data ─────────────────────────────────────────────────────────────────
tasks_download() { local r; for r in $(runs); do task download_sra "${r}" "${FASTQ_DIR}"; done; }
tasks_star() { local r; for r in $(runs); do task star_align "${r}" "$(fq "${r}" 1)" "$(fq "${r}" 2)" "${DATA_DIR}/STAR_hg38"; done; }
tasks_kallisto() { local r; for r in $(runs); do task kallisto_quant "${r}" "$(fq "${r}" 1)" "$(fq "${r}" 2)" "${DATA_DIR}/hg38_kallisto"; done; }
tasks_featurecounts() { local r; for r in $(runs); do task featurecounts_utrs "${r}" "$(bam "${r}")" "${DATA_DIR}/myUTRs_hg38"; done; }
tasks_aggregate() { task aggregate_featurecounts "${DATA_DIR}/myUTRs_hg38" df_counts.csv; }
# sum the 3 runs of each sample -> myUTRs_hg38/df_counts_summed.csv, hg38_kallisto/counts_matrix_summed.csv
tasks_sum_utr_counts() { task run_script ENV_PY --cwd "${DATA_DIR}" python "${SCRIPTS_DIR}/counts/sum_UTR_counts_across_replicates.py"; }
tasks_sum_gene_counts() { task run_script ENV_PY --cwd "${DATA_DIR}" python "${SCRIPTS_DIR}/counts/sum_genes_counts_across_replicates.py"; }

# ── genes and 3' UTRs ────────────────────────────────────────────────────────
tasks_genes_preprocessing() {
    local cfg="${RESULTS_DIR}/genes_preprocessing/genes_preprocessing.yaml"
    render_config "${DATASET_DIR}/configs/genes_preprocessing.yaml" "${cfg}"; task genes_preprocessing "${cfg}"
}
tasks_bulkapa() {
    local cfg="${RESULTS_DIR}/3UTR_analysis/bulkapa.yaml"
    render_config "${DATASET_DIR}/configs/bulkapa.yaml" "${cfg}"; task bulkapa "${cfg}"
}

# ── CLIP ─────────────────────────────────────────────────────────────────────
CLIP_RES="${RESULTS_DIR}/3UTR_analysis/CLIP"
tasks_clip_integration() {
    local kind file iso
    render_config "${DATASET_DIR}/configs/clip_eCLIP.yaml" "${CLIP_RES}/config_eCLIP.yaml"
    render_config "${DATASET_DIR}/configs/clip_iCLIP.yaml" "${CLIP_RES}/config_iCLIP.yaml"
    while read -r kind file; do
        for iso in proximal distal; do task clip_integration "${CLIP_RES}/config_${kind}.yaml" "${iso}" "${file}"; done
    done < <(list_clip_files)
}

# stop codon: same configs as clip_integration, but a -100:+100 window around the stop codon in 25-nt regions
# (results in CLIP/ next to the PAS ones, as *_results_stop_codon.pkl; drawn as the first column of the Welch profiles)
tasks_clip_stop_codon_integration() {
    local dir="${RESULTS_DIR}/3UTR_analysis/CLIP" kind file
    for kind in eCLIP iCLIP; do
        render_config "${DATASET_DIR}/configs/clip_${kind}.yaml" "${dir}/config_stop_codon_${kind}.yaml"
        override_config "${dir}/config_stop_codon_${kind}.yaml" nt_before=-100 nt_after=100 length_region=25
    done
    while read -r kind file; do task clip_integration "${dir}/config_stop_codon_${kind}.yaml" stop_codon "${file}"; done < <(list_clip_files)
}
# Welch profiles across time only (proximal, distal, combined with the stop codon) -> CLIP/analysis/; remove
# --welch-only to also redo the by-region Welch profiles and the Fisher outputs
tasks_clip_analysis() {
    render_config "${DATASET_DIR}/configs/clip_analysis.yaml" "${CLIP_RES}/clip_analysis_config.yaml"
    task clip_analysis "${CLIP_RES}/clip_analysis_config.yaml" --layout "${CLIP_ANALYSIS_LAYOUT}" --welch-only
}
tasks_clip_merge_welch() { py "${SCRIPTS_DIR}/clip/merge_welch.py" "${CLIP_RES}"; }   # -> CLIP/analysis/welch_*.csv

# Same integration with the PUD of each sample (line x day x fraction) instead of each condition (Welch tests only)
CLIP_PP_RES="${RESULTS_DIR}/3UTR_analysis/CLIP_per_patient"
tasks_clip_per_patient_integration() {
    local kind file iso
    render_config "${DATASET_DIR}/configs/clip_per_patient_eCLIP.yaml" "${CLIP_PP_RES}/config_eCLIP.yaml"
    render_config "${DATASET_DIR}/configs/clip_per_patient_iCLIP.yaml" "${CLIP_PP_RES}/config_iCLIP.yaml"
    while read -r kind file; do
        for iso in proximal distal; do task clip_integration "${CLIP_PP_RES}/config_${kind}.yaml" "${iso}" "${file}"; done
    done < <(list_clip_files)
}
# -> CLIP_per_patient/analysis/welch_*.csv
tasks_clip_per_patient_merge_welch() { py "${SCRIPTS_DIR}/clip/merge_welch_per_patient.py" "${CLIP_PP_RES}"; }

# number of array tasks, as counted by the script itself (its file list defines what --task-id points to)
script_count() { in_env "${ENV_PY}" python "$1" count "${@:2}" 2>/dev/null | tail -n 1; }

# CLIP cross-link density along whole lncRNA / paraspeckle genes, 200 nt windows
# -> paraspeckles/<GENE>/<GENE>_cross_links.csv (read by Figure_regulators_lengthening.ipynb)
PARASPECKLE_GENES=(MALAT1 NEAT1 NORAD SFPQ XIST)
para_dir() { echo "${RESULTS_DIR}/paraspeckles/$1/$1_clip_density"; }
tasks_paraspeckles_prepare() {
    local g; for g in "${PARASPECKLE_GENES[@]}"; do
        py "${SCRIPTS_DIR}/clip/gene_site_density.py" prepare --gene-name "${g}" --annotation "${GENCODE_GTF}" --output-dir "$(para_dir "${g}")"; done
}
tasks_paraspeckles_run() {
    local g i n; n=$(script_count "${SCRIPTS_DIR}/clip/gene_site_density.py" --iclip-dir "${CLIP_DIR}/iCLIP" --eclip-dir "${CLIP_DIR}/eCLIP")
    [[ "${n}" =~ ^[0-9]+$ ]] || die "could not count the CLIP files (got '${n}')"
    for g in "${PARASPECKLE_GENES[@]}"; do
        for ((i = 0; i < n; i++)); do
            py "${SCRIPTS_DIR}/clip/gene_site_density.py" run --gene-name "${g}" --iclip-dir "${CLIP_DIR}/iCLIP" --eclip-dir "${CLIP_DIR}/eCLIP" \
                --output-dir "$(para_dir "${g}")" --window-size 200 --annotation "${GENCODE_GTF}" --task-id "${i}"
        done
    done
}
tasks_paraspeckles_merge() {
    local g; for g in "${PARASPECKLE_GENES[@]}"; do
        py "${SCRIPTS_DIR}/clip/gene_site_density.py" merge --gene-name "${g}" --output-dir "$(para_dir "${g}")" --window-size 200; done
}

# ── RNA modifications ────────────────────────────────────────────────────────
RNA_MODS=(m6A m1A m5C m7G Nm Pseudo RNA-editing)
tasks_rnamod_integration() {
    local m iso cfg
    for m in "${RNA_MODS[@]}"; do
        cfg="${RESULTS_DIR}/RNA_modifs/${m}/welch/config.yaml"
        render_config "${DATASET_DIR}/configs/rnamod_${m}_welch.yaml" "${cfg}"
        for iso in distal proximal; do task clip_integration "${cfg}" "${iso}" "${RMBASE_DIR}/human.hg38.${m}.gff"; done
    done
}
# stop codon: same configs as rnamod_integration, -100:+100 window around the stop codon in 25-nt regions
# (results next to the PAS ones, as *_results_stop_codon.pkl)
tasks_rnamod_stop_codon_integration() {
    local m cfg
    for m in "${RNA_MODS[@]}"; do
        cfg="${RESULTS_DIR}/RNA_modifs/${m}/welch/config_stop_codon.yaml"
        render_config "${DATASET_DIR}/configs/rnamod_${m}_welch.yaml" "${cfg}"
        override_config "${cfg}" nt_before=-100 nt_after=100 length_region=25
        task clip_integration "${cfg}" stop_codon "${RMBASE_DIR}/human.hg38.${m}.gff"
    done
}
# Welch profiles and Fisher heatmaps of each modification -> RNA_modifs/<mod>/welch/analysis/
tasks_rnamod_analysis() {
    local m cfg
    for m in "${RNA_MODS[@]}"; do
        cfg="${RESULTS_DIR}/RNA_modifs/${m}/welch/clip_analysis_config.yaml"
        render_config "${DATASET_DIR}/configs/rnamod_${m}_analysis.yaml" "${cfg}"
        task clip_analysis "${cfg}" --layout "${CLIP_ANALYSIS_LAYOUT}"
    done
}
# -> RNA_modifs/<mod>/welch/analysis/welch_*.csv
tasks_rnamod_merge_welch() { local m; for m in "${RNA_MODS[@]}"; do py "${SCRIPTS_DIR}/clip/merge_welch.py" "${RESULTS_DIR}/RNA_modifs/${m}/welch"; done; }

# m6A structural-switch pipeline (scripts/m6a_switch/run_m6a_switch_pipeline.sh), one run per region
REGIONS_ROOT="${RESULTS_DIR}/RNA_modifs/m6A/regions"
tasks_m6a_regions() {
    local spec iso wl wr name reg
    # region 0: [-100, +100] around the 3' UTR start (stop codon), windows cut at the start of its exon; no candidate
    # list, so no CLIP test (folding and statistics only)
    task run_script ENV_BIO env PYTHON=python "EXTRA_FOLD_ARGS=--anchor utr_start --cut-at-exon-start" \
        bash "${SCRIPTS_DIR}/m6a_switch/run_m6a_switch_pipeline.sh" proximal -100 100 "${REGIONS_ROOT}/utr_start_m100_100" \
        "${M6A_SWITCH_N_JOBS}" none
    # regions 1-3; reg = region of the candidate RBP list (Tables/Figure_regulators_lengthening, from Figure_regulators_lengthening.ipynb)
    for spec in "proximal -350 0 proximal_m350_0 region1_upstream_pPAS" "proximal 0 200 proximal_0_200 region2_downstream_pPAS" "distal -350 0 distal_m350_0 region3_upstream_dPAS"; do
        read -r iso wl wr name reg <<< "${spec}"
        task run_script ENV_BIO env PYTHON=python bash "${SCRIPTS_DIR}/m6a_switch/run_m6a_switch_pipeline.sh" "${iso}" "${wl}" "${wr}" \
            "${REGIONS_ROOT}/${name}" "${M6A_SWITCH_N_JOBS}" "${REPO_ROOT}/Tables/Figure_regulators_lengthening/regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_${reg}_APA_vs_nuclear_retention_class.csv"
    done
}
tasks_m6a_regions_clip_list() { py_bio "${SCRIPTS_DIR}/m6a_switch/clip_enrichment_all_datasets.py" --root "${REGIONS_ROOT}" --list; }
tasks_m6a_regions_clip_run() {
    local k; for ((k = 0; k < 40; k++)); do
        py_bio "${SCRIPTS_DIR}/m6a_switch/clip_enrichment_all_datasets.py" --root "${REGIONS_ROOT}" --chunk "${k}" --n-chunks 40; done
}
tasks_m6a_regions_clip_merge() { py_bio "${SCRIPTS_DIR}/m6a_switch/clip_enrichment_all_datasets.py" --root "${REGIONS_ROOT}" --merge --n-chunks 40; }
tasks_m6a_regions_compare() {
    py_bio "${SCRIPTS_DIR}/m6a_switch/compare_switch_regions.py" --root "${REGIONS_ROOT}"
    py_bio "${SCRIPTS_DIR}/m6a_switch/export_pairs_switch_table.py" --root "${REGIONS_ROOT}"
}

# m6A switch competence and CLIP preference along whole lncRNAs (MALAT1, NEAT1_2, NORAD, XIST, TUG1)
LNC_DIR="${RESULTS_DIR}/RNA_modifs/m6A/lncRNA_switch"
tasks_lncrna_prepare() { py_bio "${SCRIPTS_DIR}/m6a_switch/lncrna_m6a_switch_clip.py" --out-dir "${LNC_DIR}" --prepare; }
tasks_lncrna_features() { py_bio "${SCRIPTS_DIR}/m6a_switch/lncrna_m6a_switch_clip.py" --out-dir "${LNC_DIR}" --features; }
tasks_lncrna_run() {
    local k; for ((k = 0; k < 40; k++)); do
        py_bio "${SCRIPTS_DIR}/m6a_switch/lncrna_m6a_switch_clip.py" --out-dir "${LNC_DIR}" --chunk "${k}" --n-chunks 40; done
}
tasks_lncrna_merge() { py_bio "${SCRIPTS_DIR}/m6a_switch/lncrna_m6a_switch_clip.py" --out-dir "${LNC_DIR}" --merge --n-chunks 40; }

pipeline_main "$@"
