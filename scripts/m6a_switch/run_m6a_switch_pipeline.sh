#!/usr/bin/env bash
# Full m6A switch-competence pipeline for one isoform/window, end to end:
#   1. U-tract switch-competence    (RNAfold)            -> OUT_ROOT/switch_competence_u_tract
#   2. GU-motif switch-competence   (reuses step 1 fold) -> OUT_ROOT/switch_competence_gu_motif
#   3. nuclear-lengthened vs. other statistics for steps 1 and 2 (group_stats.csv, ...)
#   4. CLIP enrichment at switch-competent U-tracts       -> OUT_ROOT/clip_u_tract_enrichment
#   5. CLIP enrichment at switch-competent GU motifs      -> OUT_ROOT/clip_gu_motif_enrichment
# Every step writes params.json (arguments, inputs, script checksums, package versions).
#
# Usage: run_m6a_switch_pipeline.sh <proximal|distal> <wl> <wr> <out_root> [n_jobs] [candidates_csv]
#   e.g. run_m6a_switch_pipeline.sh distal -350 0 /path/to/distal_m350_0 16 Tables/Figure_regulators_lengthening/regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_region3_upstream_dPAS_APA_vs_nuclear_retention_class.csv
# Without candidates_csv the CLIP steps use the scripts' default (the R2 regulators of lengthening);
# candidates_csv = none skips the CLIP steps (4-5).
# EXTRA_FOLD_ARGS is passed to step 1, e.g. EXTRA_FOLD_ARGS="--anchor utr_start --cut-at-exon-start".
# The Python interpreter can be overridden with PYTHON=/path/to/python.
set -euo pipefail

if [[ $# -lt 4 ]]; then
    sed -n '2,17p' "$0"
    exit 1
fi

ISOFORM=$1
WL=$2
WR=$3
OUT_ROOT=$4
N_JOBS=${5:-1}
CANDIDATES=${6:-}
CAND_ARGS=()
[[ -n "${CANDIDATES}" && "${CANDIDATES}" != "none" ]] && CAND_ARGS=(--candidates "$(readlink -f "${CANDIDATES}")")
read -r -a FOLD_ARGS <<< "${EXTRA_FOLD_ARGS:-}"
PYTHON=${PYTHON:-python}   # python of the bioinfo environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

U_DIR="${OUT_ROOT}/switch_competence_u_tract"
GU_DIR="${OUT_ROOT}/switch_competence_gu_motif"
CLIP_U_DIR="${OUT_ROOT}/clip_u_tract_enrichment"
CLIP_GU_DIR="${OUT_ROOT}/clip_gu_motif_enrichment"
mkdir -p "${OUT_ROOT}"

step() { echo; echo "=== [$(date '+%F %T')] $* ==="; }

step "1/5 U-tract switch-competence (${ISOFORM} [${WL}:${WR}), ${N_JOBS} jobs)"
"$PYTHON" -u "${SCRIPT_DIR}/m6a_switch_competence_analysis.py" \
    --isoform "${ISOFORM}" --wl="${WL}" --wr="${WR}" --out-dir "${U_DIR}" --n-jobs "${N_JOBS}" "${FOLD_ARGS[@]+"${FOLD_ARGS[@]}"}"

step "2/5 GU-motif switch-competence"
"$PYTHON" -u "${SCRIPT_DIR}/m6a_switch_competence_guu_motif_analysis.py" \
    --source-dir "${U_DIR}" --out-dir "${GU_DIR}"

step "3/5 group statistics"
"$PYTHON" -u "${SCRIPT_DIR}/switch_competence_group_stats.py" --run-dir "${U_DIR}"
"$PYTHON" -u "${SCRIPT_DIR}/switch_competence_group_stats.py" --run-dir "${GU_DIR}"

if [[ "${CANDIDATES}" == "none" ]]; then
    step "4-5/5 CLIP enrichment skipped (no candidate list)"
else
    step "4/5 CLIP enrichment at switch-competent U-tracts"
    "$PYTHON" -u "${SCRIPT_DIR}/clip_u_tract_switch_competent_enrichment.py" \
        --switch-dir "${U_DIR}" --out-dir "${CLIP_U_DIR}" "${CAND_ARGS[@]+"${CAND_ARGS[@]}"}"

    step "5/5 CLIP enrichment at switch-competent GU motifs"
    "$PYTHON" -u "${SCRIPT_DIR}/clip_guu_motif_switch_competent_enrichment.py" \
        --switch-dir "${GU_DIR}" --out-dir "${CLIP_GU_DIR}" "${CAND_ARGS[@]+"${CAND_ARGS[@]}"}"
fi

step "done -> ${OUT_ROOT}"
