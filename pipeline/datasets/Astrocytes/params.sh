#!/bin/bash
# Astrocytes: iPSC-derived astrocytes, VCP-mutant vs control, nuclear and cytoplasmic fractions
# (SRA PRJNA672040, 12 paired-end runs, 8 of them analysed). Samples: samples.tsv.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib.sh"

export DATASET="Astrocytes"
: "${DATASET_FOLDER:=${DATASET}}"; export DATASET_FOLDER   # output folder name in PROJECT_DIR (override to write elsewhere)
export DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/data"
export RESULTS_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/results"
export FASTQ_DIR="${DATABASE_DIR}/PRJNA672040"                # one sub-folder per run: <run>/<run>_{1,2}.fastq.gz

export STAR_INDEX="${ANNOT_DIR}/STAR/genome_index_99"
export KALLISTO_INDEX="${ANNOT_DIR}/kallisto/hg38_kallisto_v13.idx"
export KALLISTO_BOOTSTRAPS=100
export FC_STRAND=""                                           # inferred per sample with RSeQC
export STRAND_REPORT_DIR="${DATA_DIR}/rseqc_strandedness"
