#!/usr/bin/env bash
# =============================================================================
# Phase 4: Run HYMET classification on ZymoGut D6331 assembly
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
THREADS="${THREADS:-8}"
CAND_MAX="${CAND_MAX:-5000}"
SPECIES_DEDUP="${SPECIES_DEDUP:-1}"
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
Usage: 04_classify.sh [options]

Options:
  --work-dir DIR    Working directory
  --threads N       Number of threads
  --cand-max N      Maximum Mash candidates (default: 5000)
  --force           Force re-classification
  -h, --help        Show this message
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --cand-max) CAND_MAX="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

ASSEMBLY_FASTA="${WORK_DIR}/assembly/assembly.fasta"
HYMET_OUTPUT="${WORK_DIR}/hymet_output"
CLASSIFIED="${HYMET_OUTPUT}/classified_sequences.tsv"

log "ZymoGut D6331 HYMET Classification"
log "==================================="
log "Threads: ${THREADS}, Max candidates: ${CAND_MAX}"

[[ -s "${ASSEMBLY_FASTA}" ]] || die "Assembly not found: ${ASSEMBLY_FASTA}"

HYMET_CLI="${HYMET_ROOT}/bin/hymet"
[[ -x "${HYMET_CLI}" ]] || die "HYMET CLI not found: ${HYMET_CLI}"
[[ -s "${HYMET_ROOT}/data/sketch1.msh" ]] || die "HYMET sketches missing. Run: bin/hymet init"

CONTIG_COUNT=$(grep -c "^>" "${ASSEMBLY_FASTA}")
log "Input assembly: ${CONTIG_COUNT} contigs"

ELAPSED=0
if [[ -s "${CLASSIFIED}" && ${FORCE} -eq 0 ]]; then
  log "HYMET output already exists: ${CLASSIFIED}"
else
  mkdir -p "${HYMET_OUTPUT}"
  log "Running HYMET classification..."

  START_TIME=$(date +%s)
  env THREADS="${THREADS}" CAND_MAX="${CAND_MAX}" SPECIES_DEDUP="${SPECIES_DEDUP}" \
    "${HYMET_CLI}" run \
      --contigs "${ASSEMBLY_FASTA}" \
      --out "${HYMET_OUTPUT}" \
      --threads "${THREADS}" \
      --cand-max "${CAND_MAX}" \
    2>&1 | tee "${HYMET_OUTPUT}/hymet_run.log"
  END_TIME=$(date +%s)
  ELAPSED=$((END_TIME - START_TIME))
  log "HYMET completed in $((ELAPSED / 60)) minutes"
fi

[[ -s "${CLASSIFIED}" ]] || die "Classification failed"

CLASSIFIED_COUNT=$(tail -n +2 "${CLASSIFIED}" | wc -l)
UNKNOWN_COUNT=$(tail -n +2 "${CLASSIFIED}" | grep -c "Unknown" || true)
KNOWN_COUNT=$((CLASSIFIED_COUNT - UNKNOWN_COUNT))

log ""
log "Classification Results:"
log "  Total: ${CLASSIFIED_COUNT}"
log "  Classified: ${KNOWN_COUNT}"
log "  Unknown: ${UNKNOWN_COUNT}"

# Copy CAMI profile
CAMI_PROFILE="${HYMET_OUTPUT}/hymet.sample_0.cami.tsv"
PROFILE_DST="${HYMET_OUTPUT}/profile.cami.tsv"
[[ -s "${CAMI_PROFILE}" ]] && cp -f "${CAMI_PROFILE}" "${PROFILE_DST}"

log ""
log "Classification complete!"
log "Next step: Run 05_analyse.sh"

cat > "${WORK_DIR}/status_04_classify.json" <<EOF
{
  "phase": "04_classify",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "classified_sequences": "${CLASSIFIED}",
  "input_contigs": ${CONTIG_COUNT},
  "classified_count": ${CLASSIFIED_COUNT},
  "known_count": ${KNOWN_COUNT},
  "runtime_seconds": ${ELAPSED}
}
EOF
