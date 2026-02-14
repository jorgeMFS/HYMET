#!/usr/bin/env bash
# =============================================================================
# ZymoGut D6331 Analysis Workflow - Master Orchestrator
# =============================================================================
# Complete workflow for ZymoBIOMICS Gut Microbiome Standard (D6331) analysis.
# Addresses Reviewer #1's comment about ground-truth validation.
#
# Location: case/zymogut/ (self-contained case study)
#
# Phases:
#   1. Download    - Download nanopore reads and reference genomes
#   2. Assemble    - Metagenomic assembly with Flye
#   3. Build Truth - Create ground truth files
#   4. Classify    - Run HYMET classification
#   5. Analyse     - Compare predictions vs ground truth
#   6. Package     - Package results for paper integration
#
# Usage:
#   ./run_zymogut.sh                    # Run all phases
#   ./run_zymogut.sh --phase 1          # Run only phase 1
#   ./run_zymogut.sh --from 3           # Run from phase 3 onwards
#   ./run_zymogut.sh --phase 6 --publish # Package and publish to repo
# =============================================================================

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

# Defaults
WORK_DIR="${SCRIPT_DIR}/work"
THREADS="${THREADS:-8}"
CAND_MAX="${CAND_MAX:-5000}"
PHASES=""
FROM_PHASE=""
TO_PHASE=""
PUBLISH=0
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
ZymoGut D6331 Analysis Workflow
===============================

Usage: run_zymogut.sh [options]

Phase Selection:
  --phase N[,M,...]   Run specific phase(s) only (1-6)
  --from N            Start from phase N
  --to N              Stop after phase N
  (default: run all phases 1-6)

Options:
  --work-dir DIR      Working directory (default: case/zymogut/work)
  --threads N         Number of threads (default: 8 or $THREADS)
  --cand-max N        Max Mash candidates for HYMET (default: 5000)
  --publish           Publish results to repo in phase 6
  --force             Force re-run of completed phases
  -h, --help          Show this message

Phases:
  1. Download    - Download reads and references (~10 min, ~5GB)
  2. Assemble    - Flye metagenomic assembly (~30-60 min)
  3. Build Truth - Create ground truth files (~5 min)
  4. Classify    - Run HYMET classification (~10-30 min)
  5. Analyse     - Compare with ground truth (~2 min)
  6. Package     - Package results for paper (~1 min)

Examples:
  ./run_zymogut.sh                      # Full pipeline
  ./run_zymogut.sh --phase 1            # Download only
  ./run_zymogut.sh --from 4             # Re-run classification onwards
  ./run_zymogut.sh --phase 5,6 --publish # Re-analyse and publish
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --phase) PHASES="$2"; shift 2;;
    --from) FROM_PHASE="$2"; shift 2;;
    --to) TO_PHASE="$2"; shift 2;;
    --work-dir) WORK_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --cand-max) CAND_MAX="$2"; shift 2;;
    --publish) PUBLISH=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

# Determine which phases to run
declare -a RUN_PHASES=()

if [[ -n "${PHASES}" ]]; then
  IFS=',' read -ra RUN_PHASES <<< "${PHASES}"
elif [[ -n "${FROM_PHASE}" || -n "${TO_PHASE}" ]]; then
  start=${FROM_PHASE:-1}
  end=${TO_PHASE:-6}
  for ((i=start; i<=end; i++)); do
    RUN_PHASES+=("$i")
  done
else
  RUN_PHASES=(1 2 3 4 5 6)
fi

# Export common variables
export WORK_DIR THREADS CAND_MAX HYMET_ROOT CASE_ROOT SCRIPT_DIR

mkdir -p "${WORK_DIR}"

# =============================================================================
# Header
# =============================================================================
log "╔══════════════════════════════════════════════════════════════════╗"
log "║     ZymoBIOMICS Gut Microbiome Standard (D6331) Analysis         ║"
log "║                  Ground-Truth Validation Study                    ║"
log "╚══════════════════════════════════════════════════════════════════╝"
log ""
log "Configuration:"
log "  HYMET root:     ${HYMET_ROOT}"
log "  Case root:      ${CASE_ROOT}"
log "  Work directory: ${WORK_DIR}"
log "  Threads:        ${THREADS}"
log "  Phases to run:  ${RUN_PHASES[*]}"
log "  Publish:        ${PUBLISH}"
log ""

# Check HYMET is ready
if [[ ! -s "${HYMET_ROOT}/data/sketch1.msh" ]]; then
  log "WARNING: HYMET sketches not found. Run: bin/hymet init"
  log "         Or: tools/fetch_sketches.sh"
fi

START_TOTAL=$(date +%s)

# =============================================================================
# Phase execution
# =============================================================================

should_run_phase() {
  local phase=$1
  for p in "${RUN_PHASES[@]}"; do
    [[ "$p" == "$phase" ]] && return 0
  done
  return 1
}

# Phase 1: Download
if should_run_phase 1; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 1: Download"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}")
  [[ ${FORCE} -eq 1 ]] && PHASE_ARGS+=(--force)

  bash "${SCRIPT_DIR}/scripts/01_download.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# Phase 2: Assemble
if should_run_phase 2; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 2: Assembly"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}" --threads "${THREADS}")
  [[ ${FORCE} -eq 1 ]] && PHASE_ARGS+=(--force)

  bash "${SCRIPT_DIR}/scripts/02_assemble.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# Phase 3: Build Truth
if should_run_phase 3; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 3: Build Ground Truth"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}" --threads "${THREADS}")
  [[ ${FORCE} -eq 1 ]] && PHASE_ARGS+=(--force)

  bash "${SCRIPT_DIR}/scripts/03_build_truth.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# Phase 4: Classify
if should_run_phase 4; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 4: HYMET Classification"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}" --threads "${THREADS}" --cand-max "${CAND_MAX}")
  [[ ${FORCE} -eq 1 ]] && PHASE_ARGS+=(--force)

  bash "${SCRIPT_DIR}/scripts/04_classify.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# Phase 5: Analyse
if should_run_phase 5; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 5: Analysis"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}")
  [[ ${FORCE} -eq 1 ]] && PHASE_ARGS+=(--force)

  bash "${SCRIPT_DIR}/scripts/05_analyse.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# Phase 6: Package
if should_run_phase 6; then
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  log "Phase 6: Package Results"
  log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  PHASE_ARGS=(--work-dir "${WORK_DIR}")
  [[ ${PUBLISH} -eq 1 ]] && PHASE_ARGS+=(--publish)

  bash "${SCRIPT_DIR}/scripts/06_package.sh" "${PHASE_ARGS[@]}"
  log ""
fi

# =============================================================================
# Summary
# =============================================================================
END_TOTAL=$(date +%s)
ELAPSED_TOTAL=$((END_TOTAL - START_TOTAL))
ELAPSED_MIN=$((ELAPSED_TOTAL / 60))

log "╔══════════════════════════════════════════════════════════════════╗"
log "║                    Workflow Complete                              ║"
log "╚══════════════════════════════════════════════════════════════════╝"
log ""
log "Total runtime: ${ELAPSED_MIN} minutes (${ELAPSED_TOTAL} seconds)"
log ""
log "Outputs:"
log "  Work directory: ${WORK_DIR}"
log ""

# Show key results if analysis completed
if [[ -f "${WORK_DIR}/results/profile_metrics.tsv" ]]; then
  log "Key metrics:"
  while IFS=$'\t' read -r metric value; do
    [[ "$metric" == "metric" ]] && continue
    log "  ${metric}: ${value}"
  done < "${WORK_DIR}/results/profile_metrics.tsv"
  log ""
fi

if [[ -d "${WORK_DIR}/results/figures" ]]; then
  log "Figures generated:"
  for fig in "${WORK_DIR}/results/figures"/*.png; do
    [[ -f "$fig" ]] && log "  $(basename "$fig")"
  done
fi

log ""
log "Next steps:"
if [[ ${PUBLISH} -eq 0 ]]; then
  log "  1. Review results in ${WORK_DIR}/results/"
  log "  2. Run with --publish to integrate into repo:"
  log "     ./run_zymogut.sh --phase 6 --publish"
else
  log "  1. Results published to:"
  log "     - case/truth/zymogut_d6331/"
  log "     - results/cases/zymogut/run_*/"
  log "  2. Add figures to paper supplementary material"
  log "  3. Update manuscript with metrics"
fi
