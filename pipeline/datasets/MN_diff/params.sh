#!/bin/bash
# MN_diff: iPSC -> motor neuron differentiation (D0-D35), 4 control and 4 VCP-mutant lines, nuclear and cytoplasmic
# fractions (GEO GSE152983 / SRA PRJNA641118; 96 samples x 3 sequencing runs = 288 paired-end runs, summed per
# sample). Samples: samples.tsv.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib.sh"

export DATASET="MN_diff"
: "${DATASET_FOLDER:=${DATASET}}"; export DATASET_FOLDER   # output folder name in PROJECT_DIR (override to write elsewhere)
export DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/data"
export RESULTS_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/results"
export FASTQ_DIR="${DATABASE_DIR}/PRJNA641118"                # <run>_{1,2}.fastq.gz

export STAR_INDEX="${ANNOT_DIR}/STAR/genome_index_99"
export KALLISTO_INDEX="${ANNOT_DIR}/kallisto/hg38_kallisto_v13.idx"
export KALLISTO_BOOTSTRAPS=100
export FC_STRAND=""                                           # inferred per run with RSeQC
export STRAND_REPORT_DIR="${DATA_DIR}/rseqc_strandedness"

export CLIP_ANALYSIS_LAYOUT="time_course"   # Welch profiles of scripts/clip/clip_analysis.py: per day, genotype and fraction
export GENCODE_GTF="${ANNOT_DIR}/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
export M6A_SWITCH_N_JOBS=16                                   # parallel jobs inside each m6A-switch region run
