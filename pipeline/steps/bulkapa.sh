#!/bin/bash
# 3' UTR / APA analysis (apatools/scripts/bulkapa.py).
# Outputs go to the config's path_to_results (table_of_pairs.csv, shifts/, genes_utrs/, Figures/, message.log).
# Usage: bulkapa.sh CONFIG
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 1 ]] || die "usage: $(basename "$0") CONFIG"
require_file "$1"
in_env "${ENV_PY}" python "${REPOS_DIR}/apatools/scripts/bulkapa.py" --config_file "$1"
