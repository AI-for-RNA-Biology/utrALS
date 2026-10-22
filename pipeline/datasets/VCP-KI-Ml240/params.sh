#!/bin/bash
# VCP-KI-Ml240: WT / VCP-mutant / VCP knock-in / TDP-43 lines, untreated or treated with the VCP inhibitor ML240,
# nuclear and cytoplasmic fractions (SRA PRJNA883641, 47 paired-end runs). Samples: samples.tsv.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib.sh"

export DATASET="VCP-KI-Ml240"
: "${DATASET_FOLDER:=${DATASET}}"; export DATASET_FOLDER   # output folder name in PROJECT_DIR (override to write elsewhere)
export DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/data"
export RESULTS_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/results"
export FASTQ_DIR="${DATABASE_DIR}/PRJNA883641"                # one sub-folder per run: <run>/<run>_{1,2}.fastq.gz

export STAR_INDEX="${ANNOT_DIR}/STAR/genome_index_99"
export KALLISTO_INDEX="${ANNOT_DIR}/kallisto/hg38_kallisto_v13.idx"
export KALLISTO_BOOTSTRAPS=100
export FC_STRAND=""                                           # inferred per sample with RSeQC
export STRAND_REPORT_DIR="${DATA_DIR}/rseqc_strandedness"

export CLIP_ANALYSIS_LAYOUT="treatment_grid"   # Welch profiles of scripts/clip/clip_analysis.py: one line per genotype x drug
export M6A_GFF="${RMBASE_DIR}/human.hg38.m6A.gff"
