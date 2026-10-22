#!/bin/bash
# Gene-level filtering + quantile normalisation (bulkanalysis/scripts/genes_preprocessing.py).
# Outputs go to the config's path_to_results (df_genes_filtered_normalized.csv, DGE/, Figures/, message.log).
# Usage: genes_preprocessing.sh CONFIG
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 1 ]] || die "usage: $(basename "$0") CONFIG"
require_file "$1"
in_env "${ENV_PY}" python "${REPOS_DIR}/bulkanalysis/scripts/genes_preprocessing.py" --config_file "$1"
