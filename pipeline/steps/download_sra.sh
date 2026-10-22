#!/bin/bash
# Download one paired-end SRA run as gzipped FASTQ files: <FASTQ_DIR>/<RUN>_1.fastq.gz and _2.fastq.gz
# Usage: download_sra.sh RUN FASTQ_DIR
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 2 ]] || die "usage: $(basename "$0") RUN FASTQ_DIR"
RUN="$1"; FASTQ_DIR="$2"

skip_if_exists "${FASTQ_DIR}/${RUN}_2.fastq.gz"
mkdir -p "${FASTQ_DIR}"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/sra_${RUN}.XXXXXX")
trap 'rm -rf "${TMP}"' EXIT

in_env "${ENV_SRA}" fasterq-dump "${RUN}" --split-files --threads "${THREADS}" \
    --outdir "${FASTQ_DIR}" --temp "${TMP}" --skip-technical --force
require_file "${FASTQ_DIR}/${RUN}_1.fastq" "${FASTQ_DIR}/${RUN}_2.fastq"
run gzip "${FASTQ_DIR}/${RUN}_1.fastq" "${FASTQ_DIR}/${RUN}_2.fastq"
run gunzip -t "${FASTQ_DIR}/${RUN}_1.fastq.gz" "${FASTQ_DIR}/${RUN}_2.fastq.gz"
log "downloaded ${RUN}"
