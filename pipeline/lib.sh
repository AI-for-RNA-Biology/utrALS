#!/bin/bash
# Shared helpers, sourced by every pipeline script (steps, dataset runners, SLURM wrappers).

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${PIPELINE_DIR}")"
STEPS_DIR="${PIPELINE_DIR}/steps"
SCRIPTS_DIR="${REPO_ROOT}/scripts"          # command-line scripts run by the steps
# the utrals package (src/) is importable by every script, whether or not it is pip-installed in the environment
PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export PIPELINE_DIR REPO_ROOT STEPS_DIR SCRIPTS_DIR PYTHONPATH

# shellcheck source=config.sh
source "${PIPELINE_DIR}/config.sh"

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

require_file() { for f in "$@"; do [[ -e "$f" ]] || die "missing input: $f"; done; }

# Exit the calling step successfully if its output already exists (set FORCE=1 to recompute).
skip_if_exists() {
    if [[ -e "$1" && "${FORCE:-0}" != "1" ]]; then
        log "output exists, skipping (FORCE=1 to recompute): $1"
        exit 0
    fi
}

# Make the 'conda' command available (cluster module system if needed).
ensure_conda() {
    command -v conda >/dev/null 2>&1 && return 0
    if type module >/dev/null 2>&1; then module load Anaconda3 >/dev/null 2>&1 || true; fi
    command -v conda >/dev/null 2>&1 || die "conda not found; load it, or set the ENV_* variables to \"\" to use tools from PATH"
}

# in_env ENV CMD [ARGS...]: run CMD inside conda environment ENV (a name or a path).
# With ENV="" the command is run from the current PATH.
# The environment's lib/ is put first in LD_LIBRARY_PATH: otherwise the system libstdc++ (too old on the cluster)
# can be loaded instead of the environment's one (e.g. "GLIBCXX_3.4.30 not found" from llvmlite/numba).
in_env() {
    local env="$1"; shift
    log "+ ${env:+[${env}] }$*"
    if [[ -z "${env}" ]]; then
        "$@"
    else
        ensure_conda
        local flag="-n"; [[ "${env}" == */* ]] && flag="-p"
        conda run --no-capture-output "${flag}" "${env}" bash -c \
            'export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"; exec "$@"' _ "$@"
    fi
}

# run CMD [ARGS...]: log and run a command outside conda.
run() { log "+ $*"; "$@"; }

# infer_strandedness SAMPLE BAM REPORT: print the featureCounts strand flag (1 = forward, 2 = reverse, 0 = unstranded)
# from RSeQC infer_experiment.py. Rule: dominant fraction >= 0.55 and >= 2x the other one -> stranded;
# else dominant >= 0.60 -> stranded; else unstranded. An existing REPORT is reused (FORCE=1 to recompute).
infer_strandedness() {
    local sample="$1" bam="$2" report="$3" sense antisense flag
    : "${STRAND_BED:=${ANNOT_DIR}/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.bed}"
    if [[ ! -s "${report}" || "${FORCE:-0}" == "1" ]]; then
        require_file "${bam}" "${STRAND_BED}"
        mkdir -p "$(dirname "${report}")"
        in_env "${ENV_BIO}" infer_experiment.py -r "${STRAND_BED}" -i "${bam}" > "${report}"
    fi
    sense=$(grep '"1++,1--,2+-,2-+"' "${report}" | awk -F': ' '{print $2}')
    antisense=$(grep '"1+-,1-+,2++,2--"' "${report}" | awk -F': ' '{print $2}')
    flag=$(awk -v s="${sense}" -v a="${antisense}" 'BEGIN {
        if (s > a) { dom = s; oth = a; flag = 1 } else { dom = a; oth = s; flag = 2 }
        if ((dom >= 0.55 && dom >= 2 * oth) || dom >= 0.60) print flag; else print 0 }')
    log "strandedness of ${sample}: sense=${sense} antisense=${antisense} -> ${flag} (1 forward, 2 reverse, 0 unstranded)"
    echo "${flag}"
}

# Variables that may appear as ${VAR} in configs/*.yaml templates.
TEMPLATE_VARS='${PROJECT_DIR} ${DATABASE_DIR} ${ANNOT_DIR} ${CLIP_DIR} ${RMBASE_DIR} ${DATA_DIR} ${RESULTS_DIR}'

# render_config TEMPLATE OUT: fill the ${VAR} placeholders of a YAML template and write it next to the results.
render_config() {
    local template="$1" out="$2"
    require_file "${template}"
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        log "(dry-run) would render ${template} -> ${out}"
        return 0
    fi
    mkdir -p "$(dirname "${out}")"
    envsubst "${TEMPLATE_VARS}" < "${template}" > "${out}"
    if grep -n '\${' "${out}" >&2; then die "unresolved placeholder(s) in ${out}"; fi
    log "rendered ${out}"
}

# override_config FILE KEY=VALUE...: replace top-level scalar keys of a rendered config (each key must exist).
override_config() {
    local out="$1" kv key; shift
    [[ "${DRY_RUN:-0}" == "1" ]] && { log "(dry-run) would set $* in ${out}"; return 0; }
    for kv in "$@"; do
        key="${kv%%=*}"
        grep -q "^${key}:" "${out}" || die "no top-level key '${key}' in ${out}"
        sed -i "s|^${key}:.*|${key}: ${kv#*=}|" "${out}"
    done
    log "set $* in ${out}"
}

# task SCRIPT [ARGS...]: print one task line (a shell-quoted call of a step script) for a dataset runner.
task() {
    local script="$1"; shift
    printf '%q ' bash "${STEPS_DIR}/${script}.sh" "$@"
    echo
}

# list_clip_files: print "<kind> <file>" for every CLIP dataset used in the project
# (kind = eCLIP or iCLIP, which selects the config used to parse the RBP name from the file name).
# Selection identical to the original SLURM arrays: 140 eCLIP HepG2 + 177 eCLIP K562 + 18 iCLIP = 335 files.
list_clip_files() {
    local f
    for f in "${CLIP_DIR}"/eCLIP/HepG2/*.sum.gff "${CLIP_DIR}"/eCLIP/K652/*.sum.gff; do
        [[ -f "$f" ]] && echo "eCLIP $f"
    done
    local exclude="${PIPELINE_DIR}/iclip_excluded.txt"
    find "${CLIP_DIR}/iCLIP" -maxdepth 1 -type f -name "*.gff" | sort | while read -r f; do
        grep -qxF "$(basename "$f")" "${exclude}" || echo "iCLIP $f"
    done
}

# Record the code versions used for a run (utrALS + library repositories).
write_provenance() {
    local out="$1" cmdline="$2" name commit dirty
    mkdir -p "$(dirname "${out}")"
    {
        echo "date: $(date '+%F %T')"
        echo "host: $(hostname)"
        echo "command: ${cmdline}"
        dirty=$(git -C "${REPO_ROOT}" status --porcelain 2>/dev/null | wc -l)
        echo "utrALS: $(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null) (uncommitted files: ${dirty})"
        while IFS=$'\t' read -r name _url _pinned; do
            [[ "${name}" == "name" || -z "${name}" ]] && continue
            commit=$(git -C "${REPOS_DIR}/${name}" rev-parse HEAD 2>/dev/null || echo "not found")
            dirty=$(git -C "${REPOS_DIR}/${name}" status --porcelain 2>/dev/null | wc -l)
            echo "${name}: ${commit} (uncommitted files: ${dirty})"
        done < "${PIPELINE_DIR}/external_repos.tsv"
        for v in PROJECT_DIR DATABASE_DIR ANNOT_DIR CLIP_DIR RMBASE_DIR REPOS_DIR ENV_PY ENV_BIO ENV_COV ENV_SRA; do
            echo "${v}=${!v}"
        done
    } > "${out}"
}
