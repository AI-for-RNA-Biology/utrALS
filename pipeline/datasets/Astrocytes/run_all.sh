#!/bin/bash
# Regenerate PROJECT_DIR/Astrocytes (data/ and results/) from the raw FASTQ files.
#
# Usage: bash run_all.sh [--from STEP] [--to STEP] [--only STEP,...] [--jobs N] [--dry-run]
#        bash run_all.sh --list-steps
# On SLURM: bash pipeline/slurm/submit.sh Astrocytes [same step options]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/params.sh"
source "${PIPELINE_DIR}/runner.sh"

STEPS=(download star kallisto featurecounts aggregate genes_preprocessing bulkapa)

# samples: print "sample run" for every row of samples.tsv (all 12 runs are processed; the configs select 8)
samples() { awk -F'\t' 'NR > 1 && NF {print $1, $2}' "${DATASET_DIR}/samples.tsv"; }
fq() { echo "${FASTQ_DIR}/$1/$1_$2.fastq.gz"; }

tasks_download() {
    local sample run
    while read -r sample run; do task download_sra "${run}" "${FASTQ_DIR}/${run}"; done < <(samples)
}

tasks_star() {
    local sample run
    while read -r sample run; do
        task star_align "${sample}" "$(fq "${run}" 1)" "$(fq "${run}" 2)" "${DATA_DIR}/STAR_hg38"
    done < <(samples)
}

tasks_kallisto() {
    local sample run
    while read -r sample run; do
        task kallisto_quant "${sample}" "$(fq "${run}" 1)" "$(fq "${run}" 2)" "${DATA_DIR}/hg38_kallisto"
    done < <(samples)
}

tasks_featurecounts() {
    local sample run
    while read -r sample run; do
        task featurecounts_utrs "${sample}" "${DATA_DIR}/STAR_hg38/${sample}Aligned.sortedByCoord.out.bam" "${DATA_DIR}/myUTRs_hg38"
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

pipeline_main "$@"
