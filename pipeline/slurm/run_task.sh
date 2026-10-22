#!/bin/bash
# Run one line of a tasks file as a SLURM array task (called by submit.sh, not directly).
# Usage: sbatch --array=1-N run_task.sh DATASET TASKS_FILE
set -euo pipefail
DATASET="$1"; TASKS_FILE="$2"
source "${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}/datasets/${DATASET}/params.sh"

export THREADS="${SLURM_CPUS_PER_TASK:-${THREADS}}"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS_FILE}")
[[ -n "${line}" ]] || die "no task ${SLURM_ARRAY_TASK_ID} in ${TASKS_FILE}"
log "task ${SLURM_ARRAY_TASK_ID} on $(hostname): ${line}"
bash -c "${line}"
