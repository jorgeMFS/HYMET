#!/usr/bin/env bash
# =============================================================================
# Phase 2: Assemble nanopore reads with Flye
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
THREADS="${THREADS:-8}"
ASM_COVERAGE="${ASM_COVERAGE:-40}"
GENOME_SIZE="${GENOME_SIZE:-100m}"
RESUME=0
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
Usage: 02_assemble.sh [options]

Options:
  --work-dir DIR    Working directory (default: case/zymogut/work)
  --threads N       Number of threads (default: 8)
  --asm-coverage N  Target coverage for subsampling (default: 40)
  --genome-size S   Estimated metagenome size (default: 100m)
  --resume          Resume interrupted Flye assembly
  --force           Force re-assembly
  -h, --help        Show this message
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --asm-coverage) ASM_COVERAGE="$2"; shift 2;;
    --genome-size) GENOME_SIZE="$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

READS_FILE="${WORK_DIR}/downloads/zymogut_sup.fastq.gz"
ASSEMBLY_DIR="${WORK_DIR}/assembly"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly.fasta"

log "ZymoGut D6331 Assembly Script"
log "=============================="
log "Threads: ${THREADS}"
log "Assembly coverage: ${ASM_COVERAGE}x"
log "Genome size estimate: ${GENOME_SIZE}"

[[ -s "${READS_FILE}" ]] || die "Reads file not found: ${READS_FILE}\nRun 01_download.sh first."
log "Input reads: ${READS_FILE}"

# Check Flye
FLYE_BIN=""
if [[ -x "/root/.local/share/mamba/envs/hymet-benchmark/bin/flye" ]]; then
  FLYE_BIN="/root/.local/share/mamba/envs/hymet-benchmark/bin/flye"
elif command -v flye &>/dev/null; then
  FLYE_BIN="flye"
elif [[ -x "/opt/envs/squeezemeta/SqueezeMeta/bin/Flye-2.9/bin/flye" ]]; then
  FLYE_BIN="/opt/envs/squeezemeta/SqueezeMeta/bin/Flye-2.9/bin/flye"
else
  die "Flye not found. Install with: mamba install -c bioconda flye"
fi
log "Flye version: $("${FLYE_BIN}" --version 2>&1 | head -1)"

# Convert genome size to bases for subsampling calculation
parse_size() {
  local raw="$1"
  case "${raw}" in
    *[gG]) echo "$(( ${raw%[gG]} * 1000000000 ))";;
    *[mM]) echo "$(( ${raw%[mM]} * 1000000 ))";;
    *[kK]) echo "$(( ${raw%[kK]} * 1000 ))";;
    *) echo "${raw}";;
  esac
}

ELAPSED=0
if [[ -s "${ASSEMBLY_FASTA}" && ${FORCE} -eq 0 && ${RESUME} -eq 0 ]]; then
  log "Assembly already exists: ${ASSEMBLY_FASTA}"
else
  mkdir -p "${ASSEMBLY_DIR}"

  # Flye --meta does not support --asm-coverage, so subsample reads beforehand
  # to keep memory usage within bounds (~40x of estimated genome size).
  GENOME_BASES=$(parse_size "${GENOME_SIZE}")
  TARGET_BASES=$(( GENOME_BASES * ASM_COVERAGE ))
  SUBSAMPLED="${ASSEMBLY_DIR}/subsampled_reads.fastq"

  if [[ ! -s "${SUBSAMPLED}" || ${FORCE} -eq 1 ]]; then
    command -v seqtk &>/dev/null || die "seqtk not found (needed for subsampling)"
    log "Subsampling reads to ~${ASM_COVERAGE}x (~${TARGET_BASES} bp) with seqtk..."
    seqtk sample -s 42 "${READS_FILE}" "$(
      python3 -c "
import subprocess, math
# Count total bases in the input reads
result = subprocess.run(['seqtk', 'size', '${READS_FILE}'], capture_output=True, text=True)
total_bases = int(result.stdout.split()[1])
frac = min(1.0, ${TARGET_BASES} / total_bases)
print(f'{frac:.6f}')
")" > "${SUBSAMPLED}"
    SUBSAMPLE_BASES=$(seqtk size "${SUBSAMPLED}" | awk '{print $2}')
    log "Subsampled: ${SUBSAMPLE_BASES} bp"
  else
    log "Subsampled reads already exist: ${SUBSAMPLED}"
  fi

  FLYE_ARGS=(--nano-hq "${SUBSAMPLED}" --out-dir "${ASSEMBLY_DIR}" --meta --threads "${THREADS}")
  [[ ${RESUME} -eq 1 ]] && FLYE_ARGS+=(--resume)

  # Ensure minimap2 is on PATH (Flye needs it internally)
  MINIMAP2_DIR="$(dirname "$(command -v minimap2 2>/dev/null || echo /opt/conda/bin/minimap2)")"
  export PATH="${MINIMAP2_DIR}:${PATH}"

  log "Running Flye (this may take 30-60 minutes)..."
  START_TIME=$(date +%s)
  "${FLYE_BIN}" "${FLYE_ARGS[@]}" 2>&1 | tee "${ASSEMBLY_DIR}/flye.log"
  END_TIME=$(date +%s)
  ELAPSED=$((END_TIME - START_TIME))
  log "Assembly completed in $((ELAPSED / 60)) minutes"
fi

[[ -s "${ASSEMBLY_FASTA}" ]] || die "Assembly failed"

CONTIG_COUNT=$(grep -c "^>" "${ASSEMBLY_FASTA}")
TOTAL_LENGTH=$(awk '/^>/{next}{total+=length}END{print total}' "${ASSEMBLY_FASTA}")

log ""
log "Assembly Statistics:"
log "  Contigs: ${CONTIG_COUNT}"
log "  Total length: ${TOTAL_LENGTH} bp"

log ""
log "Assembly complete!"
log "Next step: Run 03_build_truth.sh"

cat > "${WORK_DIR}/status_02_assemble.json" <<EOF
{
  "phase": "02_assemble",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "assembly_fasta": "${ASSEMBLY_FASTA}",
  "contig_count": ${CONTIG_COUNT},
  "total_length": ${TOTAL_LENGTH},
  "runtime_seconds": ${ELAPSED}
}
EOF
