#!/bin/bash
# Align one paired-end sample with STAR -> <OUTDIR>/<SAMPLE>Aligned.sortedByCoord.out.bam
# Usage: STAR_INDEX=<genome dir> star_align.sh SAMPLE R1 R2 OUTDIR
#   R1/R2 can be comma-separated lists (several lanes of the same sample).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 4 ]] || die "usage: STAR_INDEX=<dir> $(basename "$0") SAMPLE R1 R2 OUTDIR"
SAMPLE="$1"; R1="$2"; R2="$3"; OUTDIR="$4"
: "${STAR_INDEX:?set STAR_INDEX (STAR genome directory) in the dataset params.sh}"

skip_if_exists "${OUTDIR}/${SAMPLE}Aligned.sortedByCoord.out.bam"
require_file ${R1//,/ } ${R2//,/ } "${STAR_INDEX}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"
STDERR="./${SAMPLE}STAR.stderr"
run "${STAR_BIN}" --runThreadN "${THREADS}" --genomeDir "${STAR_INDEX}" --readFilesIn "${R1}" "${R2}" \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "./${SAMPLE}" --readFilesCommand zcat \
    2> >(tee "${STDERR}" >&2)

# STAR can exit successfully while warning that the final BAM is corrupt ("Unexpected block structure ... Possible
# output corruption" while merging the sorted bins). Seen for SFPQ_TDP_FUS_KD sample TYZ4001A14, reproducibly
# (2 runs, 2026-10-08/09), although the original BAM of the same sample is fine. Remove such a BAM and fail.
BAM="${SAMPLE}Aligned.sortedByCoord.out.bam"
if grep -q "Possible output corruption" "${STDERR}"; then
    rm -f "${BAM}"
    die "STAR reported possible output corruption for ${SAMPLE}; BAM removed, rerun this task"
fi
in_env "${ENV_COV}" samtools quickcheck "${BAM}" || { rm -f "${BAM}"; die "truncated BAM for ${SAMPLE}; removed"; }
rm -f "${STDERR}"
