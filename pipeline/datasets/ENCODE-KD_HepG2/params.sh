#!/bin/bash
# ENCODE-KD_HepG2: ENCODE shRNA knockdowns of single RBPs in HepG2 cells, each compared with the shared control
# (SRR4422171, SRR4422172). Runs per group (control or RBP): samples.tsv. Steps: ../encode_kd.sh.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib.sh"

export DATASET="ENCODE-KD_HepG2"
: "${DATASET_FOLDER:=${DATASET}}"; export DATASET_FOLDER   # output folder name in PROJECT_DIR (override to write elsewhere)
export DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export RESULTS_DIR="${PROJECT_DIR}/${DATASET_FOLDER}"                # one sub-folder per RBP
# Raw processing (FASTQ, STAR, kallisto, featureCounts) is kept outside the project folder:
: "${ENCODE_DATA_DIR:=${EXTERNAL_DATA_DIR}/HepG2_shRNA}"
export ENCODE_DATA_DIR
export DATA_DIR="${ENCODE_DATA_DIR}"

export STAR_INDEX="${ANNOT_DIR}/STAR/genome_index_99"
export KALLISTO_INDEX="${ANNOT_DIR}/kallisto/hg38_kallisto.idx"
export KALLISTO_BOOTSTRAPS=0
export KALLISTO_STRAND="auto"                                 # inferred per run with RSeQC
export FC_STRAND=""                                           # inferred per run with RSeQC
