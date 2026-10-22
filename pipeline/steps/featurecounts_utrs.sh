#!/bin/bash
# Count reads on the 3' UTR isoform annotation for one sample -> <OUTDIR>/myUTRs_hg38_<SAMPLE>.txt
# Usage: [FC_STRAND=0|1|2] [UTR_GTF=<gtf>] featurecounts_utrs.sh SAMPLE BAM OUTDIR
#   If FC_STRAND is empty, the strandedness is inferred with RSeQC infer_experiment.py
#   (report in ${STRAND_REPORT_DIR:-<OUTDIR>/../rseqc_strandedness}/<SAMPLE>_infer_experiment.txt):
#     dominant fraction >= 0.55 and >= 2x the other one -> stranded; else >= 0.60 -> stranded; else unstranded.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 3 ]] || die "usage: $(basename "$0") SAMPLE BAM OUTDIR"
SAMPLE="$1"; BAM="$2"; OUTDIR="$3"
: "${UTR_GTF:=${ANNOT_DIR}/3UTR_MNs/hg38_NoOverlap.gtf}"
: "${FC_STRAND:=}"

OUT="${OUTDIR}/myUTRs_hg38_${SAMPLE}.txt"
skip_if_exists "${OUT}"
require_file "${BAM}" "${UTR_GTF}"
mkdir -p "${OUTDIR}"

if [[ -z "${FC_STRAND}" ]]; then
    FC_STRAND=$(infer_strandedness "${SAMPLE}" "${BAM}" \
        "${STRAND_REPORT_DIR:-$(dirname "${OUTDIR}")/rseqc_strandedness}/${SAMPLE}_infer_experiment.txt")
fi

in_env "${ENV_BIO}" featureCounts -a "${UTR_GTF}" -o "${OUT}" -g transcript_id -f -O -M \
    -s "${FC_STRAND}" -p -T "${THREADS}" --verbose "${BAM}"
