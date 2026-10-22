#!/bin/bash
# Run a project script (Python, R, shell) inside one of the pipeline's conda environments, optionally from a given
# working directory. Used for the dataset-specific analyses that have no dedicated step.
# Usage: run_script.sh ENV_VARIABLE [--cwd DIR] COMMAND [ARGS...]
#   ENV_VARIABLE: name of the config.sh variable holding the environment (ENV_PY, ENV_BIO, ...)
#   e.g. run_script.sh ENV_PY python scripts/clip/merge_welch.py <folder>
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -ge 2 ]] || die "usage: $(basename "$0") ENV_VARIABLE [--cwd DIR] COMMAND [ARGS...]"
ENV_VAR="$1"; shift
[[ -v "${ENV_VAR}" ]] || die "unknown environment variable ${ENV_VAR} (see pipeline/config.sh)"
if [[ "$1" == "--cwd" ]]; then
    mkdir -p "$2"; cd "$2"; shift 2
fi
in_env "${!ENV_VAR}" "$@"
