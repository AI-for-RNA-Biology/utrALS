#!/bin/bash
# Generic dataset runner, sourced at the end of each pipeline/datasets/<DATASET>/run_all.sh.
#
# The dataset's run_all.sh defines:
#   STEPS=(step1 step2 ...)              the ordered list of steps
#   tasks_<step>() { task <script> args; ... }
#                                        prints one line per independent task of that step
#                                        (e.g. one per sample, one per CLIP file); may render configs.
# and then calls: pipeline_main "$@"
#
# Options:
#   --from STEP / --to STEP   run a range of steps (default: all)
#   --only STEP[,STEP...]     run only these steps
#   --jobs N                  run up to N tasks of a step in parallel (default 1)
#   --dry-run                 print the tasks without running them
#   --list-tasks STEP         print the task lines of one step (used by pipeline/slurm/submit.sh)
#   --list-steps              print the step names
#   --print-steps             print the steps selected by --from/--to/--only (used by pipeline/slurm/submit.sh)

usage() { sed -n '/^# Options:/,/^$/p' "${PIPELINE_DIR}/runner.sh" | sed 's/^# \{0,1\}//'; }

# select_steps FROM TO ONLY: print the selected steps, in pipeline order.
select_steps() {
    local from="$1" to="$2" only="$3" s on=0
    if [[ -n "${only}" ]]; then
        for s in "${STEPS[@]}"; do
            [[ ",${only}," == *",${s},"* ]] && echo "${s}"
        done
        return
    fi
    [[ -z "${from}" ]] && on=1
    for s in "${STEPS[@]}"; do
        [[ "${s}" == "${from}" ]] && on=1
        (( on )) && echo "${s}"
        [[ "${s}" == "${to}" ]] && break
    done
}

check_step() {
    local s
    for s in "${STEPS[@]}"; do [[ "${s}" == "$1" ]] && return 0; done
    die "unknown step '$1'. Steps: ${STEPS[*]}"
}

# run_step STEP: generate the step's tasks and run them locally, N at a time, one log file per task.
run_step() {
    local step="$1" tasks_file n i=0 failed=0 line
    tasks_file=$(mktemp)
    "tasks_${step}" > "${tasks_file}" || die "could not generate the tasks of step ${step}"
    n=$(grep -c . "${tasks_file}" || true)
    log "== step ${step}: ${n} task(s)"
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        cat "${tasks_file}"; rm -f "${tasks_file}"; return 0
    fi
    local log_dir="${LOG_DIR}/${step}"
    mkdir -p "${log_dir}"
    while IFS= read -r line; do
        [[ -z "${line}" ]] && continue
        i=$((i + 1))
        while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || failed=$((failed + 1)); done
        (
            if bash -c "${line}" > "${log_dir}/task_${i}.log" 2>&1; then
                log "   ${step} ${i}/${n} done"
            else
                log "   ${step} ${i}/${n} FAILED (see ${log_dir}/task_${i}.log)"; exit 1
            fi
        ) &
    done < "${tasks_file}"
    while (( $(jobs -rp | wc -l) > 0 )); do wait -n || failed=$((failed + 1)); done
    rm -f "${tasks_file}"
    (( failed == 0 )) || die "step ${step}: ${failed} task(s) failed"
}

pipeline_main() {
    local from="" to="" only="" list_tasks="" s print_steps=0
    local cmdline="$0 $*"
    JOBS=1
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --from) from="$2"; check_step "$2"; shift 2 ;;
            --to) to="$2"; check_step "$2"; shift 2 ;;
            --only) only="$2"; for s in ${2//,/ }; do check_step "$s"; done; shift 2 ;;
            --jobs) JOBS="$2"; shift 2 ;;
            --dry-run) export DRY_RUN=1; shift ;;
            --list-tasks) list_tasks="$2"; check_step "$2"; shift 2 ;;
            --list-steps) printf '%s\n' "${STEPS[@]}"; return 0 ;;
            --print-steps) print_steps=1; shift ;;
            -h|--help) usage; return 0 ;;
            *) usage; die "unknown option $1" ;;
        esac
    done

    if [[ -n "${list_tasks}" ]]; then
        "tasks_${list_tasks}"
        return
    fi
    if (( print_steps )); then
        select_steps "${from}" "${to}" "${only}"
        return
    fi

    LOG_DIR="${PROJECT_DIR}/${DATASET_FOLDER}/pipeline_logs/$(date +%Y%m%d_%H%M%S)"
    if [[ "${DRY_RUN:-0}" != "1" ]]; then
        write_provenance "${LOG_DIR}/provenance.txt" "${cmdline}"
        log "logs and code versions: ${LOG_DIR}"
    fi
    for s in $(select_steps "${from}" "${to}" "${only}"); do
        run_step "${s}"
    done
    log "finished"
}
