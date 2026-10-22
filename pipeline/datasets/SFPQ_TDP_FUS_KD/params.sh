#!/bin/bash
# SFPQ_TDP_FUS_KD: siRNA knockdown of SFPQ, TDP-43 and FUS (+ untreated, scramble), nuclear and cytoplasmic
# fractions, 3 replicates (30 samples, in-house run RN21152, paired-end over 2 lanes). Samples: samples.tsv.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib.sh"

export DATASET="SFPQ_TDP_FUS_KD"
: "${DATASET_FOLDER:=${DATASET}}"; export DATASET_FOLDER   # output folder name in PROJECT_DIR (override to write elsewhere)
export DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/data"
export RESULTS_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/results"
export FASTQ_DIR="${DATABASE_DIR}/SFPQ_dataset/RN21152/210924_A01366_0060_AH3CYFDMXY/fastq"

export STAR_INDEX="${ANNOT_DIR}/STAR/genome_index_99"
export KALLISTO_INDEX="${ANNOT_DIR}/kallisto/hg38_kallisto_v13.idx"
export KALLISTO_BOOTSTRAPS=100
export FC_STRAND=""                                          # inferred per sample with RSeQC
export STRAND_REPORT_DIR="${DATA_DIR}/rseqc_strandedness"

export CLIP_ANALYSIS_LAYOUT="knockdown"   # Welch profiles of scripts/clip/clip_analysis.py: one line per treatment
export M6A_GFF="${RMBASE_DIR}/human.hg38.m6A.gff"
