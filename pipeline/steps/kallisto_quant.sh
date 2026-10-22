#!/bin/bash
# Quantify one paired-end sample with kallisto -> <OUTDIR>/<SAMPLE>/abundance.tsv
# Usage: KALLISTO_INDEX=<idx> [KALLISTO_BOOTSTRAPS=100] [KALLISTO_STRAND=--rf-stranded] kallisto_quant.sh SAMPLE R1 R2 OUTDIR [BAM]
#   R1/R2 can be comma-separated lists (several lanes): they are passed as R1_a R2_a R1_b R2_b ... to kallisto.
#   KALLISTO_STRAND=auto: infer the strandedness from BAM with RSeQC (report in
#   ${STRAND_REPORT_DIR:-<OUTDIR>/../rseqc_strandedness}/<SAMPLE>_infer_experiment.txt) -> --fr-stranded / --rf-stranded / none.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"

[[ $# -eq 4 || $# -eq 5 ]] || die "usage: KALLISTO_INDEX=<idx> $(basename "$0") SAMPLE R1 R2 OUTDIR [BAM]"
SAMPLE="$1"; R1="$2"; R2="$3"; OUTDIR="$4"; BAM="${5:-}"
: "${KALLISTO_INDEX:?set KALLISTO_INDEX in the dataset params.sh}"
: "${KALLISTO_BOOTSTRAPS:=0}"
: "${KALLISTO_STRAND=--rf-stranded}"

skip_if_exists "${OUTDIR}/${SAMPLE}/abundance.tsv"
IFS=',' read -r -a r1 <<< "${R1}"
IFS=',' read -r -a r2 <<< "${R2}"
[[ ${#r1[@]} -eq ${#r2[@]} ]] || die "R1 and R2 have a different number of files"
reads=()
for i in "${!r1[@]}"; do reads+=("${r1[$i]}" "${r2[$i]}"); done
require_file "${reads[@]}" "${KALLISTO_INDEX}"

if [[ "${KALLISTO_STRAND}" == "auto" ]]; then
    [[ -n "${BAM}" ]] || die "KALLISTO_STRAND=auto needs the BAM file as 5th argument"
    case $(infer_strandedness "${SAMPLE}" "${BAM}" "${STRAND_REPORT_DIR:-$(dirname "${OUTDIR}")/rseqc_strandedness}/${SAMPLE}_infer_experiment.txt") in
        1) KALLISTO_STRAND="--fr-stranded" ;;
        2) KALLISTO_STRAND="--rf-stranded" ;;
        *) KALLISTO_STRAND="" ;;
    esac
fi

opts=(--index "${KALLISTO_INDEX}" -o "${OUTDIR}/${SAMPLE}" -t "${THREADS}")
[[ -n "${KALLISTO_STRAND}" ]] && opts+=("${KALLISTO_STRAND}")
(( KALLISTO_BOOTSTRAPS > 0 )) && opts+=(-b "${KALLISTO_BOOTSTRAPS}")
mkdir -p "${OUTDIR}"
run "${KALLISTO_BIN}" quant "${opts[@]}" "${reads[@]}"
