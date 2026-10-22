#!/bin/bash
# Regenerate PROJECT_DIR/SFPQ_TDP_FUS_KD (data/ and results/) from the raw FASTQ files.
#
# Usage: bash run_all.sh [--from STEP] [--to STEP] [--only STEP,...] [--jobs N] [--dry-run]
#        bash run_all.sh --list-steps
# On SLURM: bash pipeline/slurm/submit.sh SFPQ_TDP_FUS_KD [same step options]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/params.sh"
source "${PIPELINE_DIR}/runner.sh"

STEPS=(star kallisto featurecounts aggregate genes_preprocessing bulkapa
       clip_integration clip_stop_codon_integration clip_analysis m6a_integration m6a_stop_codon_integration m6a_analysis)

# samples: print "sample condition fastq_prefix" for every row of samples.tsv
samples() { awk -F'\t' 'NR > 1 && NF {print $1, $2, $3}' "${DATASET_DIR}/samples.tsv"; }

# Two lanes per sample: R1 = L001_R1,L002_R1 and R2 = L001_R2,L002_R2 (paired lane by lane)
r1() { echo "${FASTQ_DIR}/$1_L001_R1_001.fastq.gz,${FASTQ_DIR}/$1_L002_R1_001.fastq.gz"; }
r2() { echo "${FASTQ_DIR}/$1_L001_R2_001.fastq.gz,${FASTQ_DIR}/$1_L002_R2_001.fastq.gz"; }
bam() { echo "${DATA_DIR}/STAR_hg38/$1Aligned.sortedByCoord.out.bam"; }

tasks_star() {
    local sample condition prefix
    while read -r sample condition prefix; do
        task star_align "${sample}" "$(r1 "${prefix}")" "$(r2 "${prefix}")" "${DATA_DIR}/STAR_hg38"
    done < <(samples)
}

tasks_kallisto() {
    local sample condition prefix
    while read -r sample condition prefix; do
        task kallisto_quant "${sample}" "$(r1 "${prefix}")" "$(r2 "${prefix}")" "${DATA_DIR}/hg38_kallisto"
    done < <(samples)
}

tasks_featurecounts() {
    local sample condition prefix
    while read -r sample condition prefix; do
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

# one analysis per knockdown -> CLIP/analysis_{SFPQ,TDP43,FUS}
tasks_clip_analysis() {
    local kd cfg
    for kd in SFPQ TDP43 FUS; do
        cfg="${RESULTS_DIR}/3UTR_analysis/CLIP/clip_analysis_config_${kd}.yaml"
        render_config "${DATASET_DIR}/configs/clip_analysis_${kd}.yaml" "${cfg}"
        task clip_analysis "${cfg}" --layout "${CLIP_ANALYSIS_LAYOUT}"
    done
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

# one analysis per knockdown -> RNA_modifs/m6A/welch/analysis_{SFPQ,TDP43,FUS}
tasks_m6a_analysis() {
    local kd cfg
    for kd in SFPQ TDP43 FUS; do
        cfg="${RESULTS_DIR}/RNA_modifs/m6A/welch/clip_analysis_${kd}.yaml"
        render_config "${DATASET_DIR}/configs/m6a_analysis_${kd}.yaml" "${cfg}"
        task clip_analysis "${cfg}" --layout "${CLIP_ANALYSIS_LAYOUT}"
    done
}

pipeline_main "$@"
