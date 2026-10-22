#!/bin/bash
# Regenerate the ENCODE K562 knockdown analyses (PROJECT_DIR/ENCODE-KD_K562/<RBP>/) from the raw FASTQ files.
#
# Usage: bash run_all.sh [--from STEP] [--to STEP] [--only STEP,...] [--jobs N] [--dry-run]
#        ENCODE_RBPS="TARDBP FUS" bash run_all.sh ...      (only some RBPs; the control is always included)
# On SLURM: bash pipeline/slurm/submit.sh ENCODE-KD_K562 [same step options]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/params.sh"
source "${PIPELINE_DIR}/runner.sh"
source "${PIPELINE_DIR}/datasets/encode_kd.sh"

pipeline_main "$@"
