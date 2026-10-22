#!/bin/bash
# Steps shared by the ENCODE shRNA knockdown datasets (ENCODE-KD_K562, ENCODE-KD_HepG2), sourced by their run_all.sh.
#
# Each knocked-down RBP is compared with the cell line's shared control (2 runs). Raw processing is stored per group
# (control or RBP) under ENCODE_DATA_DIR/<group>/{raw,STAR_hg38,kallisto_hg38,myUTRs_hg38,rseqc_strandedness};
# the analyses go to PROJECT_DIR/<dataset>/<RBP>/{genes_preprocessing,3UTR_analysis}.
# Set ENCODE_RBPS="TARDBP FUS ..." to process only some RBPs (default: all RBPs of samples.tsv).
# The analysis configs are written as the original per-RBP generator scripts did (all existing configs regenerate identically).

STEPS=(download star kallisto featurecounts genes_preprocessing bulkapa)

# samples: print "run group" for the control and the selected RBPs
samples() {
    awk -F'\t' -v sel=" ${ENCODE_RBPS:-} " 'NR > 1 && NF {
        if ($2 == "control" || sel == "  " || index(sel, " " $2 " ")) print $1, $2 }' "${DATASET_DIR}/samples.tsv"
}
rbps() { samples | awk '$2 != "control" && !seen[$2]++ {print $2}'; }
runs_of() { samples | awk -v g="$1" '$2 == g {print $1}'; }
fq() { echo "${ENCODE_DATA_DIR}/$2/raw/$1_$3.fastq.gz"; }   # run group 1|2
bam() { echo "${ENCODE_DATA_DIR}/$2/STAR_hg38/$1Aligned.sortedByCoord.out.bam"; }

tasks_download() {
    local run group
    while read -r run group; do task download_sra "${run}" "${ENCODE_DATA_DIR}/${group}/raw"; done < <(samples)
}

tasks_star() {
    local run group
    while read -r run group; do
        task star_align "${run}" "$(fq "${run}" "${group}" 1)" "$(fq "${run}" "${group}" 2)" "${ENCODE_DATA_DIR}/${group}/STAR_hg38"
    done < <(samples)
}

# strandedness inferred from the BAM (KALLISTO_STRAND=auto in params.sh)
tasks_kallisto() {
    local run group
    while read -r run group; do
        task kallisto_quant "${run}" "$(fq "${run}" "${group}" 1)" "$(fq "${run}" "${group}" 2)" \
            "${ENCODE_DATA_DIR}/${group}/kallisto_hg38" "$(bam "${run}" "${group}")"
    done < <(samples)
}

tasks_featurecounts() {
    local run group
    while read -r run group; do
        task featurecounts_utrs "${run}" "$(bam "${run}" "${group}")" "${ENCODE_DATA_DIR}/${group}/myUTRs_hg38"
    done < <(samples)
}

write_genes_preprocessing_config() {
    local rbp="$1" out="${RESULTS_DIR}/$1/genes_preprocessing/genes_preprocessing.yaml" run
    if [[ "${DRY_RUN:-0}" == "1" ]]; then log "(dry-run) would write ${out}"; echo "${out}"; return; fi
    mkdir -p "$(dirname "${out}")"
    {
        echo "data_origin: kallisto_whole_transcriptome"
        echo "df_counts_path: ${RESULTS_DIR}/${rbp}/genes_preprocessing"
        echo "data_list: "
        local g; for g in control "${rbp}"; do
            for run in $(runs_of "${g}"); do echo "  - ${ENCODE_DATA_DIR}/${g}/kallisto_hg38/${run}/abundance.tsv"; done
        done
        echo "gene_names_path: \"${ANNOT_DIR}/GENCODE/Gencode_geneNames_hg38V44.txt\""
        echo "gene_info_path: \"${ANNOT_DIR}/Homo_sapiens.gene_info\""
        echo "treatments:"
        echo "  control: "
        for run in $(runs_of control); do echo "    - '${run}'"; done
        echo "  ${rbp}-KD:"
        for run in $(runs_of "${rbp}"); do echo "    - '${run}'"; done
        echo "path_to_results: \"${RESULTS_DIR}/${rbp}/genes_preprocessing\""
        echo "figures_extension: \"pdf\""
        echo "pct_in_treatment: 0.75"
        echo "remove_outliers: False"
    } > "${out}"
    log "wrote ${out}"
    echo "${out}"
}

write_bulkapa_config() {
    local rbp="$1" out="${RESULTS_DIR}/$1/3UTR_analysis/bulkapa.yaml" run
    if [[ "${DRY_RUN:-0}" == "1" ]]; then log "(dry-run) would write ${out}"; echo "${out}"; return; fi
    mkdir -p "$(dirname "${out}")"
    {
        echo "df_counts_3UTR_path: \"${RESULTS_DIR}/${rbp}/3UTR_analysis/df_counts.csv\""
        echo "data_list: "
        local g; for g in control "${rbp}"; do
            for run in $(runs_of "${g}"); do echo "  - \"${ENCODE_DATA_DIR}/${g}/myUTRs_hg38/myUTRs_hg38_${run}.txt\""; done
        done
        echo "gtf_file_focus: \"${ANNOT_DIR}/3UTR_MNs/hg38_NoOverlap.gtf\""
        echo "gtf_file_no_focus: \"${ANNOT_DIR}/3UTR_MNs/hg38_Denovo_Ensembl_RefSeq_Filt_Reformat.gtf\""
        echo "df_filtered_genes_path: \"${RESULTS_DIR}/${rbp}/genes_preprocessing/df_genes_filtered_normalized.csv\""
        echo "figures_extension: \"pdf\""
        echo "pct_in_treatment: 0.5"
        echo "quantile_threshold: 0.05"
        echo "background_quantile: False"
        echo "treatments:"
        echo "  control: "
        for run in $(runs_of control); do echo "    - '${run}'"; done
        echo "  ${rbp}-KD:"
        for run in $(runs_of "${rbp}"); do echo "    - '${run}'"; done
        echo "thresh_FDR: 0.01"
        echo "thresh_pval: 0.05"
        echo "thresh_PUD: 0.15"
        echo "thresh_RUD: 1"
        echo "path_to_results: \"${RESULTS_DIR}/${rbp}/3UTR_analysis\""
    } > "${out}"
    log "wrote ${out}"
    echo "${out}"
}

tasks_genes_preprocessing() {
    local rbp
    for rbp in $(rbps); do task genes_preprocessing "$(write_genes_preprocessing_config "${rbp}")"; done
}

tasks_bulkapa() {
    local rbp
    for rbp in $(rbps); do task bulkapa "$(write_bulkapa_config "${rbp}")"; done
}
