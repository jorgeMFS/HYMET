#!/usr/bin/env bash
# =============================================================================
# Phase 1: Download ZymoGut D6331 data
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
SKIP_READS=0
SKIP_REFS=0
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
Usage: 01_download.sh [options]

Options:
  --work-dir DIR    Working directory (default: case/zymogut/work)
  --skip-reads      Skip downloading nanopore reads
  --skip-refs       Skip downloading reference genomes
  --force           Force re-download even if files exist
  -h, --help        Show this message
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --skip-reads) SKIP_READS=1; shift;;
    --skip-refs) SKIP_REFS=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

READS_URL="http://ftp.sra.ebi.ac.uk/vol1/run/ERR142/ERR14251410/PAY12289_barcode13.dorado0.8.2.bm5.0.0.sup.sim.fastq.gz"
REFS_URL="https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip"

DOWNLOAD_DIR="${WORK_DIR}/downloads"
mkdir -p "${DOWNLOAD_DIR}"

log "ZymoGut D6331 Download Script"
log "=============================="
log "Work directory: ${WORK_DIR}"
echo

# Download nanopore reads
READS_FILE="${DOWNLOAD_DIR}/zymogut_sup.fastq.gz"

if [[ ${SKIP_READS} -eq 0 ]]; then
  if [[ -s "${READS_FILE}" && ${FORCE} -eq 0 ]]; then
    log "Reads file already exists: ${READS_FILE}"
    log "  Size: $(du -h "${READS_FILE}" | cut -f1)"
  else
    log "Downloading nanopore SUP basecalls..."
    log "  URL: ${READS_URL}"
    wget --progress=bar:force:noscroll -O "${READS_FILE}.tmp" "${READS_URL}"
    mv "${READS_FILE}.tmp" "${READS_FILE}"
    log "Downloaded: ${READS_FILE}"
  fi
  log "Reads sanity check:"
  zcat "${READS_FILE}" | head -4 || die "Cannot read FASTQ file"
  READS_COUNT=$(zcat "${READS_FILE}" | awk 'NR%4==1' | wc -l)
  log "  Read count: ${READS_COUNT}"
else
  log "Skipping reads download (--skip-reads)"
  READS_COUNT=0
fi

echo

# Download reference genomes
REFS_ZIP="${DOWNLOAD_DIR}/D6331.refseq.zip"
REFS_DIR="${DOWNLOAD_DIR}/zymogut_refs"

if [[ ${SKIP_REFS} -eq 0 ]]; then
  if [[ -d "${REFS_DIR}" && -n "$(ls -A "${REFS_DIR}" 2>/dev/null)" && ${FORCE} -eq 0 ]]; then
    log "Reference directory already exists: ${REFS_DIR}"
  else
    log "Downloading reference genomes..."
    wget --progress=bar:force:noscroll -O "${REFS_ZIP}.tmp" "${REFS_URL}"
    mv "${REFS_ZIP}.tmp" "${REFS_ZIP}"
    rm -rf "${REFS_DIR}"
    mkdir -p "${REFS_DIR}"
    unzip -q "${REFS_ZIP}" -d "${REFS_DIR}"
    if [[ -d "${REFS_DIR}/D6331.refseq" ]]; then
      mv "${REFS_DIR}/D6331.refseq"/* "${REFS_DIR}/" 2>/dev/null || true
      rmdir "${REFS_DIR}/D6331.refseq" 2>/dev/null || true
    fi
    log "Extracted to: ${REFS_DIR}"
  fi
  log "Reference genomes:"
  shopt -s nullglob
  for f in "${REFS_DIR}"/*.fasta "${REFS_DIR}"/*.fa "${REFS_DIR}"/*.fna; do
    [[ -f "$f" ]] && log "  $(basename "$f")"
  done
  shopt -u nullglob
else
  log "Skipping references download (--skip-refs)"
fi

log ""
log "Download complete!"
log "Next step: Run 02_assemble.sh"

cat > "${WORK_DIR}/status_01_download.json" <<EOF
{
  "phase": "01_download",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "reads_file": "${READS_FILE}",
  "refs_dir": "${REFS_DIR}",
  "reads_count": ${READS_COUNT:-0}
}
EOF
