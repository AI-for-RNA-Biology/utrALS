#!/bin/bash
# Submit the steps of a dataset to SLURM: one job array per step (one array task per sample / CLIP file),
# each step waiting for the previous one to succeed. The tasks are exactly those that the dataset's
# run_all.sh runs locally.
#
# Usage: bash pipeline/slurm/submit.sh DATASET [--from STEP] [--to STEP] [--only STEP,...] [--dry-run]
#   --dry-run   print the sbatch commands without submitting
# Run it with 'bash' from a login node (it calls sbatch itself). Settings: pipeline/slurm/slurm.conf
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib.sh"
source "${PIPELINE_DIR}/slurm/slurm.conf"

[[ $# -ge 1 ]] || die "usage: $(basename "$0") DATASET [--from STEP] [--to STEP] [--only STEP,...] [--dry-run]"
DATASET="$1"; shift
RUN_ALL="${PIPELINE_DIR}/datasets/${DATASET}/run_all.sh"
require_file "${RUN_ALL}"

dry=0; select_opts=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) dry=1; shift ;;
        --from|--to|--only) select_opts+=("$1" "$2"); shift 2 ;;
        *) die "unknown option $1" ;;
    esac
done

source "${PIPELINE_DIR}/datasets/${DATASET}/params.sh"
if (( dry )); then
    SUBMIT_DIR=$(mktemp -d)
else
    SUBMIT_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/pipeline_logs/slurm_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "${SUBMIT_DIR}"
fi
(( dry )) || write_provenance "${SUBMIT_DIR}/provenance.txt" "$0 ${DATASET} ${select_opts[*]}"

previous=""
for step in $(bash "${RUN_ALL}" --print-steps "${select_opts[@]}"); do
    tasks="${SUBMIT_DIR}/${step}.tasks"
    if (( dry )); then
        DRY_RUN=1 bash "${RUN_ALL}" --list-tasks "${step}" > "${tasks}"
    else
        bash "${RUN_ALL}" --list-tasks "${step}" > "${tasks}"
    fi
    n=$(grep -c . "${tasks}" || true)
    if (( n == 0 )); then log "step ${step}: no task, skipped"; continue; fi

    res_var="RES_${step}"; read -r cpus mem time <<< "${!res_var:-${RES_default}}"
    cmd=(sbatch --parsable --job-name="${DATASET}_${step}" --array="1-${n}%${SLURM_MAX_PARALLEL}"
         --nodes=1 --cpus-per-task="${cpus}" --mem="${mem}" --time="${time}" --account="${SLURM_ACCOUNT}"
         --output="${SUBMIT_DIR}/${step}_%a.out" --error="${SUBMIT_DIR}/${step}_%a.err")
    [[ -n "${previous}" ]] && cmd+=(--dependency="afterok:${previous}")
    # shellcheck disable=SC2206
    [[ -n "${SLURM_EXTRA_OPTS}" ]] && cmd+=(${SLURM_EXTRA_OPTS})
    cmd+=("${PIPELINE_DIR}/slurm/run_task.sh" "${DATASET}" "${tasks}")

    if (( dry )); then
        echo "${cmd[*]}"
        previous="<${step}_jobid>"
    else
        previous=$("${cmd[@]}")
        log "step ${step}: ${n} task(s) -> job ${previous}"
    fi
done
if (( dry )); then rm -rf "${SUBMIT_DIR}"; else log "logs: ${SUBMIT_DIR}"; fi
