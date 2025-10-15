#!/usr/bin/env bash
# Run HYMET classifier on a sample (wrapper for run_hymet_cami.sh).

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"

usage(){
  cat <<'USAGE'
Usage: run_hymet.sh --sample ID --contigs FASTA [--out DIR] [--threads N]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/hymet}"
ensure_dir "${OUT_DIR}"
RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "Input contigs FASTA missing (${CONTIGS_ABS})"
fi

log "Running HYMET classifier for ${SAMPLE}"
INPUT_FASTA="${CONTIGS_ABS}" \
OUTDIR="${RUN_DIR}" \
THREADS="${THREADS}" \
ROOT="${HYMET_ROOT}" \
bash "${HYMET_ROOT}/run_hymet_cami.sh"

PROFILE_SRC="${RUN_DIR}/hymet.sample_0.cami.tsv"
CLASSIFIED_SRC="${RUN_DIR}/classified_sequences.tsv"
PAF_SRC="${RUN_DIR}/work/resultados.paf"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
PAF_DST="${OUT_DIR}/resultados.paf"

if [[ -s "${PROFILE_SRC}" ]]; then
  cp -f "${PROFILE_SRC}" "${PROFILE_DST}"
else
  die "Expected HYMET CAMI profile not found at ${PROFILE_SRC}"
fi

CLASSIFIED_META=""
if [[ -s "${CLASSIFIED_SRC}" ]]; then
  cp -f "${CLASSIFIED_SRC}" "${CLASSIFIED_DST}"
  CLASSIFIED_META="${CLASSIFIED_DST}"
fi

PAF_META=""
if [[ -s "${PAF_SRC}" ]]; then
  cp -f "${PAF_SRC}" "${PAF_DST}"
  PAF_META="${PAF_DST}"
fi

if [[ "${KEEP_HYMET_WORK:-0}" -eq 0 ]]; then
  rm -f "${RUN_DIR}/work/reference.mmi"
  rm -f "${CLASSIFIED_SRC}" "${PAF_SRC}" "${PROFILE_SRC}"
fi

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "hymet", "profile": "${PROFILE_DST}", "contigs": "${CLASSIFIED_META}", "paf": "${PAF_META}", "input_fasta": "${CONTIGS_ABS}", "run_dir": "${RUN_DIR}"}
EOF
