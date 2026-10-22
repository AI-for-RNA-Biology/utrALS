#!/bin/bash
# Global configuration of the pipeline: paths, tools and conda environments.
# Every value can be overridden by exporting the variable before running, e.g.
#   PROJECT_DIR=/my/scratch bash pipeline/datasets/Astrocytes/run_all.sh
# The data locations (PROJECT_DIR, DATABASE_DIR, ANNOT_DIR, CLIP_DIR, RMBASE_DIR, REPOS_DIR, ...) are defined in
# paths.env at the repository root, shared with the notebooks and scripts (utrals.paths).

# ── Data locations: paths.env, an exported variable taking precedence ────────
if [[ ! -f "${REPO_ROOT}/paths.env" ]]; then
    echo "ERROR: ${REPO_ROOT}/paths.env not found: copy paths.env.example to paths.env and set DATA_ROOT." >&2
    exit 1
fi
while IFS='=' read -r _key _value || [[ -n "${_key}" ]]; do
    [[ "${_key}" =~ ^[A-Z_][A-Z0-9_]*$ ]] || continue          # comments and blank lines
    [[ -n "${!_key:-}" ]] || eval "${_key}=\"${_value}\""     # ${KEY} refers to a variable defined above
    export "${_key}"
done < "${REPO_ROOT}/paths.env"
unset _key _value
# PROJECT_DIR: where results are written (one sub-folder per dataset); DATABASE_DIR: raw FASTQ files;
# ANNOT_DIR: genome annotations / indices; CLIP_DIR: eCLIP/{HepG2,K652}, iCLIP .gff; RMBASE_DIR: RNA modification
# sites .gff; REPOS_DIR: clones of apatools and bulkanalysis (scripts run by the pipeline) at the commits listed in
# pipeline/external_repos.tsv (run pipeline/setup_external_repos.sh to clone them).

# ── Tools not installed through conda ─────────────────────────────────────────
_default_tool() { command -v "$1" 2>/dev/null || echo "$2"; }
: "${STAR_BIN:=$(_default_tool STAR "${GENOMICS_DIR}/bin/STAR/source/STAR")}"
: "${KALLISTO_BIN:=$(_default_tool kallisto "${GENOMICS_DIR}/bin/kallisto-lnx_amd64_v0.46.1/kallisto")}"

# ── Conda environments (see envs/*.yml). Set to "" to use the tools from the current PATH ──
: "${ENV_PY=utrals}"        # Python analyses: bulkapa, genes_preprocessing, CLIP integration, MN_diff scripts
: "${ENV_BIO=bioinfo}"      # featureCounts, RSeQC infer_experiment.py; MN_diff m6A-switch scripts
: "${ENV_COV=coverage_env}" # samtools (BAM check after STAR)
: "${ENV_SRA=sra_download}" # sra-tools fasterq-dump

# ── Resources ─────────────────────────────────────────────────────────────────
: "${THREADS:=8}"           # threads per task (SLURM tasks use --cpus-per-task instead)

export STAR_BIN KALLISTO_BIN ENV_PY ENV_BIO ENV_COV ENV_SRA THREADS
