#!/bin/bash
# Regenerate PROJECT_DIR/VCP-KI-Ml240 (data/ and results/) from the raw FASTQ files.
#
# Usage: bash run_all.sh [--from STEP] [--to STEP] [--only STEP,...] [--jobs N] [--dry-run]
#        bash run_all.sh --list-steps
# On SLURM: bash pipeline/slurm/submit.sh VCP-KI-Ml240 [same step options]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/params.sh"
source "${PIPELINE_DIR}/runner.sh"

STEPS=(download star kallisto featurecounts aggregate genes_preprocessing bulkapa
       clip_integration clip_stop_codon_integration clip_analysis m6a_integration m6a_stop_codon_integration m6a_analysis)

# samples: print "sample condition" for every row of samples.tsv (sample = SRA run)
samples() { awk -F'\t' 'NR > 1 && NF {print $1, $2}' "${DATASET_DIR}/samples.tsv"; }
fq() { echo "${FASTQ_DIR}/$1/$1_$2.fastq.gz"; }
bam() { echo "${DATA_DIR}/STAR_hg38/$1Aligned.sortedByCoord.out.bam"; }

tasks_download() {
    local sample condition
    while read -r sample condition; do task download_sra "${sample}" "${FASTQ_DIR}/${sample}"; done < <(samples)
}

tasks_star() {
    local sample condition
    while read -r sample condition; do
        task star_align "${sample}" "$(fq "${sample}" 1)" "$(fq "${sample}" 2)" "${DATA_DIR}/STAR_hg38"
    done < <(samples)
}

tasks_kallisto() {
    local sample condition
    while read -r sample condition; do
        task kallisto_quant "${sample}" "$(fq "${sample}" 1)" "$(fq "${sample}" 2)" "${DATA_DIR}/hg38_kallisto"
    done < <(samples)
}

tasks_featurecounts() {
    local sample condition
    while read -r sample condition; do
        task featurecounts_utrs "${sample}" "$(bam "${sample}")" "${DATA_DIR}/myUTRs_hg38"
    done < <(samples)
}

tasks_aggregate() { task aggregate_featurecounts "${DATA_DIR}/myUTRs_hg38" df_counts.csv; }

tasks_genes_preprocessing() {
    local cfg="${RESULTS_DIR}/genes_preprocessing/genes_preprocessing.yaml"
    render_config "${DATASET_DIR}/configs/genes_preprocessing.yaml" "${cfg}"
    task genes_preprocessing "${cfg}"
}

tasks_bulkapa() {
    local cfg="${RESULTS_DIR}/3UTR_analysis/bulkapa.yaml"
    render_config "${DATASET_DIR}/configs/bulkapa.yaml" "${cfg}"
    task bulkapa "${cfg}"
}

# one task per CLIP file and isoform type (335 files x 2)
tasks_clip_integration() {
    local dir="${RESULTS_DIR}/3UTR_analysis/CLIP" kind file iso
    render_config "${DATASET_DIR}/configs/clip_eCLIP.yaml" "${dir}/config_eCLIP.yaml"
    render_config "${DATASET_DIR}/configs/clip_iCLIP.yaml" "${dir}/config_iCLIP.yaml"
    while read -r kind file; do
        for iso in proximal distal; do task clip_integration "${dir}/config_${kind}.yaml" "${iso}" "${file}"; done
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

# Welch profiles (Nuclear / Cytoplasmic x stop codon / proximal / distal) -> CLIP/analysis/
tasks_clip_analysis() {
    local cfg="${RESULTS_DIR}/3UTR_analysis/CLIP/clip_analysis_config.yaml"
    render_config "${DATASET_DIR}/configs/clip_analysis.yaml" "${cfg}"
    task clip_analysis "${cfg}" --layout "${CLIP_ANALYSIS_LAYOUT}"
}

tasks_m6a_integration() {
    local cfg="${RESULTS_DIR}/RNA_modifs/m6A/welch/config.yaml" iso
    render_config "${DATASET_DIR}/configs/m6a_welch.yaml" "${cfg}"
    for iso in distal proximal; do task clip_integration "${cfg}" "${iso}" "${M6A_GFF}"; done
}

# stop codon: same config as m6a_integration, -100:+100 window around the stop codon in 25-nt regions
# (results next to the PAS ones, as *_results_stop_codon.pkl)
tasks_m6a_stop_codon_integration() {
    local cfg="${RESULTS_DIR}/RNA_modifs/m6A/welch/config_stop_codon.yaml"
    render_config "${DATASET_DIR}/configs/m6a_welch.yaml" "${cfg}"
    override_config "${cfg}" nt_before=-100 nt_after=100 length_region=25
    task clip_integration "${cfg}" stop_codon "${M6A_GFF}"
}

# m6A Welch profiles -> RNA_modifs/m6A/welch/analysis/
tasks_m6a_analysis() {
    local cfg="${RESULTS_DIR}/RNA_modifs/m6A/welch/welch_plot_config.yaml"
    render_config "${DATASET_DIR}/configs/m6a_analysis.yaml" "${cfg}"
    task clip_analysis "${cfg}" --layout "${CLIP_ANALYSIS_LAYOUT}"
}

pipeline_main "$@"
