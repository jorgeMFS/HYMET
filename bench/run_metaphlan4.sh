#!/usr/bin/env bash
# Run MetaPhlAn4 on contigs and convert output to CAMI format.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
THREADS="${METAPHLAN_THREADS:-${THREADS}}"
EXTRA_OPTS="${METAPHLAN_OPTS:-}"
TAXDIR="${TAXDIR:-${HYMET_ROOT}/taxonomy_files}"
METAPHLAN_CMD="${METAPHLAN_CMD:-metaphlan}"
DB_DIR="${METAPHLAN_DB_DIR:-${BENCH_ROOT}/db/metaphlan}"
INDEX_NAME="${METAPHLAN_INDEX:-mpa_vJun23_CHOCOPhlAnSGB_202307}"

usage(){
  cat <<'USAGE'
Usage: run_metaphlan4.sh --sample ID --contigs FASTA [--out DIR] [--threads N] [--opts "..."]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --opts) EXTRA_OPTS="$2"; shift 2;;
    --taxdir) TAXDIR="$2"; shift 2;;
    --cmd) METAPHLAN_CMD="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/metaphlan4}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "MetaPhlAn input FASTA missing (${CONTIGS_ABS})"
fi

if ! command -v "${METAPHLAN_CMD}" >/dev/null 2>&1; then
  die "${METAPHLAN_CMD} executable not found. Please install MetaPhlAn4 and set METAPHLAN_CMD if needed."
fi

if [[ ! -d "${DB_DIR}" ]]; then
  die "MetaPhlAn database directory missing (${DB_DIR})."
fi

if [[ ! -s "${DB_DIR}/${INDEX_NAME}.pkl" ]]; then
  die "MetaPhlAn index ${INDEX_NAME} not found under ${DB_DIR}."
fi

RAW_PROFILE="${OUT_DIR}/metaphlan_profile.tsv"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
IFS=' ' read -r -a EXTRA_ARGS <<< "${EXTRA_OPTS:-}"

log "Running MetaPhlAn4 for ${SAMPLE}"
CMD=(
  "${METAPHLAN_CMD}"
  "${CONTIGS_ABS}"
  "--input_type" "fasta"
  "--nproc" "${THREADS}"
  "--db_dir" "${DB_DIR}"
  "-x" "${INDEX_NAME}"
)
if [[ ${#EXTRA_ARGS[@]} -gt 0 && -n "${EXTRA_ARGS[0]}" ]]; then
  CMD+=("${EXTRA_ARGS[@]}")
fi
CMD+=("-o" "${RAW_PROFILE}")
"${CMD[@]}"

python3 "${SCRIPT_DIR}/convert/metaphlan4_to_cami.py" \
  --input "${RAW_PROFILE}" \
  --out "${PROFILE_CAMI}" \
  --sample-id "${SAMPLE}" \
  --tool "metaphlan4" \
  --taxdb "${TAXDIR}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "metaphlan4", "profile": "${PROFILE_CAMI}", "raw_profile": "${RAW_PROFILE}", "metaphlan_cmd": "${METAPHLAN_CMD}", "extra_opts": "${EXTRA_OPTS}", "db_dir": "${DB_DIR}", "index": "${INDEX_NAME}"}
EOF
