#!/usr/bin/env bash
# =============================================================================
# Phase 6: Package ZymoGut D6331 results for paper integration
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
PUBLISH=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --publish) PUBLISH=1; shift;;
    -h|--help) echo "Usage: 06_package.sh [--work-dir DIR] [--publish]"; exit 0;;
    *) shift;;
  esac
done

TRUTH_DIR="${WORK_DIR}/truth"
RESULTS_DIR="${WORK_DIR}/results"
FIGURES_DIR="${RESULTS_DIR}/figures"
HYMET_OUTPUT="${WORK_DIR}/hymet_output"
ASSEMBLY_DIR="${WORK_DIR}/assembly"

CASE_TRUTH_DIR="${CASE_ROOT}/truth/zymogut_d6331"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RESULTS_RUN_DIR="${HYMET_ROOT}/results/cases/zymogut/run_${RUN_STAMP}"

log "ZymoGut D6331 Results Packaging"
log "================================"

[[ -s "${TRUTH_DIR}/truth_contigs.tsv" ]] || die "Truth files not found"
[[ -s "${RESULTS_DIR}/comparison_table.tsv" ]] || die "Analysis results not found"

log "Files to package:"
log "  Truth -> ${CASE_TRUTH_DIR}/"
log "  Results -> ${RESULTS_RUN_DIR}/"
log ""

if [[ ${PUBLISH} -eq 0 ]]; then
  log "DRY RUN - No files copied."
  log "Use --publish to actually copy files."
  exit 0
fi

# Copy truth files
log "Copying truth files..."
mkdir -p "${CASE_TRUTH_DIR}"
cp -v "${TRUTH_DIR}/truth_contigs.tsv" "${CASE_TRUTH_DIR}/"
cp -v "${TRUTH_DIR}/truth_profile.cami.tsv" "${CASE_TRUTH_DIR}/"
cp -v "${TRUTH_DIR}/seqid2taxid.tsv" "${CASE_TRUTH_DIR}/"
cp -v "${TRUTH_DIR}/contigs_vs_refs.paf" "${CASE_TRUTH_DIR}/zymogut_d6331_vs_refs.paf"

# Create results structure
log "Creating results structure..."
RAW_DIR="${RESULTS_RUN_DIR}/raw/zymogut_d6331"
TABLES_DIR="${RESULTS_RUN_DIR}/tables"
FIGS_DIR="${RESULTS_RUN_DIR}/figures"
mkdir -p "${RAW_DIR}" "${TABLES_DIR}" "${FIGS_DIR}"

cp -v "${HYMET_OUTPUT}/classified_sequences.tsv" "${RAW_DIR}/"
cp -v "${HYMET_OUTPUT}/profile.cami.tsv" "${RAW_DIR}/" 2>/dev/null || \
  cp -v "${HYMET_OUTPUT}/hymet.sample_0.cami.tsv" "${RAW_DIR}/profile.cami.tsv"
[[ -f "${HYMET_OUTPUT}/metadata.json" ]] && cp -v "${HYMET_OUTPUT}/metadata.json" "${RAW_DIR}/"

cp -v "${RESULTS_DIR}/comparison_table.tsv" "${TABLES_DIR}/"
cp -v "${RESULTS_DIR}/profile_metrics.tsv" "${TABLES_DIR}/"
cp -v "${RESULTS_DIR}/analysis_summary.json" "${TABLES_DIR}/"
cp -v "${FIGURES_DIR}"/*.png "${FIGS_DIR}/" 2>/dev/null || true

# Create metadata
L1=$(grep "l1_distance" "${RESULTS_DIR}/profile_metrics.tsv" 2>/dev/null | cut -f2 || echo "N/A")
BC=$(grep "bray_curtis" "${RESULTS_DIR}/profile_metrics.tsv" 2>/dev/null | cut -f2 || echo "N/A")
CORR=$(grep "correlation" "${RESULTS_DIR}/profile_metrics.tsv" 2>/dev/null | cut -f2 || echo "N/A")

cat > "${RESULTS_RUN_DIR}/metadata.json" <<EOF
{
  "sample_id": "zymogut_d6331",
  "dataset": "ZymoBIOMICS Gut Microbiome Standard (D6331)",
  "run_id": "run_${RUN_STAMP}",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "source": {
    "reads": "MicroBench (Kirk3gaard) - Oxford Nanopore SUP basecalls",
    "references": "Zymo Research D6331.refseq.zip"
  },
  "metrics": {
    "l1_distance": ${L1},
    "bray_curtis": ${BC},
    "correlation": ${CORR}
  }
}
EOF

# Create manifest entry
MANIFEST_ZYMOGUT="${CASE_ROOT}/manifest_zymogut.tsv"
cat > "${MANIFEST_ZYMOGUT}" <<EOF
sample_id	contigs_fa	truth_contigs_tsv	truth_profile_tsv	expected_taxa	citation
zymogut_d6331	${ASSEMBLY_DIR}/assembly.fasta	truth/zymogut_d6331/truth_contigs.tsv	truth/zymogut_d6331/truth_profile.cami.tsv	"Faecalibacterium prausnitzii;Veillonella rogosae;Bacteroides fragilis;Roseburia hominis;Lactobacillus fermentum;Bifidobacterium adolescentis;Fusobacterium nucleatum;Saccharomyces cerevisiae;Candida albicans;Prevotella corporis;Clostridioides difficile;Akkermansia muciniphila"	"ZymoBIOMICS Gut Microbiome Standard (D6331); MicroBench (Kirk3gaard)"
EOF

log ""
log "================================"
log "Packaging complete!"
log ""
log "Published to:"
log "  Truth:    ${CASE_TRUTH_DIR}/"
log "  Results:  ${RESULTS_RUN_DIR}/"
log "  Manifest: ${MANIFEST_ZYMOGUT}"
log ""
log "Figures for paper:"
shopt -s nullglob
for fig in "${FIGS_DIR}"/*.png; do
  [[ -f "$fig" ]] && log "  - $(basename "$fig")"
done
shopt -u nullglob

cat > "${WORK_DIR}/status_06_package.json" <<EOF
{
  "phase": "06_package",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "case_truth_dir": "${CASE_TRUTH_DIR}",
  "results_run_dir": "${RESULTS_RUN_DIR}"
}
EOF
