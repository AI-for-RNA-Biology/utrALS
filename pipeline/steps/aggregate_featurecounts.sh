#!/bin/bash
# Merge the per-sample featureCounts outputs of a folder (myUTRs_hg38_<SAMPLE>.txt) into one matrix.
# Usage: aggregate_featurecounts.sh DIR [OUT_NAME=df_counts.csv]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -ge 1 ]] || die "usage: $(basename "$0") DIR [OUT_NAME]"
DIR="$1"; OUT_NAME="${2:-df_counts.csv}"

cd "${DIR}"
files=(myUTRs_hg38_*.txt)
[[ -e "${files[0]}" ]] || die "no myUTRs_hg38_*.txt file in ${DIR}"
names=()
for f in "${files[@]}"; do n="${f#myUTRs_hg38_}"; names+=("${n%.txt}"); done
in_env "${ENV_PY}" python "${REPOS_DIR}/bulkanalysis/scripts/aggregate_featureCounts_output.py" \
    -f "${files[@]}" -n "${names[@]}" -s "${OUT_NAME}"
