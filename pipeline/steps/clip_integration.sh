#!/bin/bash
# Integrate one CLIP (or RMBase modification) .gff file with the APA pairs, around the proximal or distal PAS, or around the stop codon
# (apatools/scripts/clip_integration.py) -> <path_to_results>/<RBP>_<folder>_results_<ISO_TYPE>.pkl
# Usage: clip_integration.sh CONFIG ISO_TYPE(proximal|distal|stop_codon) GFF
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 3 ]] || die "usage: $(basename "$0") CONFIG proximal|distal|stop_codon GFF"
CONFIG="$1"; ISO_TYPE="$2"; GFF="$3"
[[ "${ISO_TYPE}" =~ ^(proximal|distal|stop_codon)$ ]] || die "ISO_TYPE must be proximal, distal or stop_codon"
require_file "${CONFIG}" "${GFF}"
in_env "${ENV_PY}" python "${REPOS_DIR}/apatools/scripts/clip_integration.py" --config_file "${CONFIG}" \
    --iso_type "${ISO_TYPE}" --n_jobs 1 --clip_file "${GFF}"
