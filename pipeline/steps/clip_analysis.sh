#!/bin/bash
# Welch profiles (and Fisher plots) of all integrated CLIP / RNA-modification results of a folder, with
# scripts/clip/clip_analysis.py; the layout depends on the dataset (CLIP_ANALYSIS_LAYOUT in its params.sh).
# Usage: clip_analysis.sh CONFIG --layout {time_course,knockdown,treatment_grid} [--welch-only]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -ge 3 ]] || die "usage: $(basename "$0") CONFIG --layout LAYOUT [--welch-only]"
require_file "$1"
in_env "${ENV_PY}" python "${SCRIPTS_DIR}/clip/clip_analysis.py" "$@"
