#!/bin/bash
# Clone the library repositories used by the pipeline into REPOS_DIR, at the commits pinned in
# external_repos.tsv. Repositories that already exist are never modified: if their checked-out
# code differs from the pinned commit, a warning is printed (later commits that only update poetry.lock are accepted).
#
# Usage: bash pipeline/setup_external_repos.sh            (REPOS_DIR from config.sh or the environment)
#        REPOS_DIR=/my/repos bash pipeline/setup_external_repos.sh

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

mkdir -p "${REPOS_DIR}"
status=0
while IFS=$'\t' read -r name url commit; do
    [[ "${name}" == "name" || -z "${name}" ]] && continue
    dest="${REPOS_DIR}/${name}"
    if [[ -d "${dest}/.git" ]]; then
        current=$(git -C "${dest}" rev-parse HEAD)
        if [[ "${current}" == "${commit}" ]]; then
            log "${name}: OK (${commit:0:7})"
        elif git -C "${dest}" diff --quiet "${commit}" HEAD -- . ':(exclude)poetry.lock' 2>/dev/null; then
            log "${name}: OK (${current:0:7}: pinned ${commit:0:7} + poetry.lock updates)"
        else
            log "WARNING ${name}: checked out at ${current:0:7}, pinned commit is ${commit:0:7} (not changed)"
            status=1
        fi
        if [[ -n "$(git -C "${dest}" status --porcelain)" ]]; then
            log "WARNING ${name}: has uncommitted changes"
            status=1
        fi
    else
        log "${name}: cloning ${url} at ${commit:0:7}"
        git clone --quiet "${url}" "${dest}"
        git -C "${dest}" checkout --quiet "${commit}"
    fi
done < "${PIPELINE_DIR}/external_repos.tsv"

exit "${status}"
