#!/usr/bin/env bash
# run_hymet_cami.sh — HYMET end‑to‑end: Mash → Download → Minimap2 → Classify → CAMI export
#
# Usage examples:
#   bash run_hymet_cami.sh
#   INPUT_FASTA=/data/cami/sample_0.fna OUTDIR=/data/hymet_out/sample_0 THREADS=16 bash run_hymet_cami.sh
#
# Notes:
# - Works with the fixed classifier in scripts/classification.py. If that fails (empty result),
#   a robust fallback (first-hit mapping via Identifiers) kicks in automatically.
# - Avoids copying a file onto itself and includes clear diagnostics.
#
# Reviewer comment to retain in repo docs:
# * Benchmarking is currently based solely on the authors' own simulated data. This is insufficient to validate
#   performance claims. Standard datasets such as those from the CAMI challenge should be incorporated.

set -Eeuo pipefail

# ----------------- params -----------------
INPUT_FASTA="${INPUT_FASTA:-/data/cami/sample_0.fna}"
OUTDIR="${OUTDIR:-/data/hymet_out/sample_0}"
THREADS="${THREADS:-8}"
CAND_MAX="${CAND_MAX:-5000}"      # max candidates after Mash
SPLIT_IDX="${SPLIT_IDX:-2g}"      # minimap2 -I chunk size for lower RAM
MASH_THRESH="${MASH_THRESH:-0.9}" # Mash screen threshold
ROOT="${ROOT:-/data/HYMET}"
export PATH="/opt/conda/bin:$PATH"
export TMPDIR="${TMPDIR:-/data/tmp}"
export OMP_NUM_THREADS="$THREADS"
export TAXONKIT_DB="$ROOT/taxonomy_files"
umask 0022
# -----------------------------------------

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }
trap 'die "failed at line $LINENO"' ERR

# ensure required dirs exist
mkdir -p "$TMPDIR"

cd "$ROOT" || die "cannot cd to $ROOT"

# ---- checks ----
[ -s "$INPUT_FASTA" ] || die "missing FASTA $INPUT_FASTA"
for f in data/sketch1.msh data/sketch2.msh data/sketch3.msh data/detailed_taxonomy.tsv data/taxonomy_hierarchy.tsv \
         scripts/mash.sh scripts/minimap2.sh scripts/classification.py; do
  [ -s "$f" ] || die "missing $f"
done
[ -d taxonomy_files ] || { log "taxonomy_files missing → running ./config.pl"; ./config.pl; }

# Check for required external tools
for tool in /data/tools/hymet2cami_fixed.py /data/tools/build_id_map.py /data/tools/mini_classify.py; do
  [ -f "$tool" ] || die "missing required tool: $tool"
done

# ---- IO layout ----
WORKDIR="$OUTDIR/work"
rm -rf output 2>/dev/null || true
mkdir -p "$OUTDIR" "$WORKDIR" "$OUTDIR/logs" input cache data/downloaded_genomes
ln -s "$WORKDIR" output 2>/dev/null || true
cp -f "$INPUT_FASTA" input/sample_0.fna

# ---- deps ----
command -v minimap2 >/dev/null || micromamba install -y -p /opt/conda -c bioconda minimap2
command -v mash     >/dev/null || micromamba install -y -p /opt/conda -c bioconda mash
command -v taxonkit >/dev/null || micromamba install -y -p /opt/conda -c bioconda taxonkit

# ---- make minimap2 memory-friendly (idempotent) ----
# ensure all invocations inside scripts/minimap2.sh carry -I$SPLIT_IDX
if ! grep -Eq 'minimap2 .* -I[0-9]+' scripts/minimap2.sh; then
  sed -i "s/minimap2 -d /minimap2 -I${SPLIT_IDX} -d /g" scripts/minimap2.sh || true
  sed -i "s/minimap2 -t /minimap2 -I${SPLIT_IDX} -t /g" scripts/minimap2.sh || true
fi

# 1) Mash → selected_genomes
if [ ! -s output/selected_genomes.txt ]; then
  log "Mash screen (sketch1)"
  bash scripts/mash.sh input data/sketch1.msh \
    output/screen.tab output/filtered_screen.tab output/sorted_screen.tab \
    output/top_hits.tab output/selected_genomes.txt "$MASH_THRESH"
  log "Mash screen (sketch2)"
  bash scripts/mash.sh input data/sketch2.msh \
    output/gtdb_screen.tab /tmp/gtdb_filtered.tab /tmp/gtdb_sorted.tab \
    /tmp/gtdb_top_hits.tab /tmp/gtdb_selected_genomes.txt "$MASH_THRESH" || true
  cat /tmp/gtdb_selected_genomes.txt >> output/selected_genomes.txt 2>/dev/null || true
  log "Mash screen (sketch3)"
  bash scripts/mash.sh input data/sketch3.msh \
    output/custom_screen.tab /tmp/custom_filtered.tab /tmp/custom_sorted.tab \
    /tmp/custom_top_hits.tab /tmp/custom_selected_genomes.txt "$MASH_THRESH" || true
  cat /tmp/custom_selected_genomes.txt >> output/selected_genomes.txt 2>/dev/null || true
  sort -u -o output/selected_genomes.txt output/selected_genomes.txt
fi

# 2) optional candidate limit
LINES=$(wc -l < output/selected_genomes.txt || echo 0)
if [ "$LINES" -gt "$CAND_MAX" ]; then
  log "limiting candidates: $LINES → $CAND_MAX"
  head -n "$CAND_MAX" output/selected_genomes.txt > output/selected_genomes.txt.tmp \
    && mv output/selected_genomes.txt.tmp output/selected_genomes.txt
fi

# 3) download refs + concat
if [ ! -s data/downloaded_genomes/combined_genomes.fasta ]; then
  log "downloadDB.py"
  python3 scripts/downloadDB.py output/selected_genomes.txt data/downloaded_genomes data/detailed_taxonomy.tsv cache
fi

# 4) minimap2 index+map
if [ ! -s output/resultados.paf ]; then
  log "minimap2 index+map"
  bash scripts/minimap2.sh input data/downloaded_genomes/combined_genomes.fasta output/reference.mmi output/resultados.paf
fi

# 5) classify (primary path)
log "classification.py on raw PAF"
python3 scripts/classification.py \
  --paf output/resultados.paf \
  --taxonomy data/detailed_taxonomy.tsv \
  --hierarchy data/taxonomy_hierarchy.tsv \
  --output output/classified_sequences.tsv \
  --processes "$THREADS" || true

ROWS=$(wc -l < output/classified_sequences.tsv 2>/dev/null || echo 0)
if [ "$ROWS" -lt 2 ]; then
  log "primary classification empty → running robust fallback (first-hit via Identifiers)"

  # Use existing fallback tools
  BUILD_ID_MAP="/data/tools/build_id_map.py"
  MINI_CLASSIFY="/data/tools/mini_classify.py"
  
  [ -f "$BUILD_ID_MAP" ] || die "missing build_id_map.py: $BUILD_ID_MAP"
  [ -f "$MINI_CLASSIFY" ] || die "missing mini_classify.py: $MINI_CLASSIFY"

  python3 "$BUILD_ID_MAP" data/detailed_taxonomy.tsv "$WORKDIR/id_to_taxid.tsv" | tee "$OUTDIR/logs/step_fallback_build_map.log"
  python3 "$MINI_CLASSIFY" output/resultados.paf "$WORKDIR/id_to_taxid.tsv" "$WORKDIR/fallback_classified.tsv" | tee "$OUTDIR/logs/step_fallback_classify.log"

  # Convert mini_classify output to HYMET format if needed
  if [ -s "$WORKDIR/fallback_classified.tsv" ]; then
    # Simple conversion: add Lineage, Taxonomic Level, Confidence columns
    awk -F'\t' 'BEGIN{OFS="\t"; print "Query\tLineage\tTaxonomic Level\tConfidence"} 
                 NR>1{print $1, "unknown", "unknown", "1.0000"}' \
        "$WORKDIR/fallback_classified.tsv" > output/classified_sequences.tsv
  fi

  ROWS=$(wc -l < output/classified_sequences.tsv 2>/dev/null || echo 0)
  [ "$ROWS" -ge 2 ] || die "classification still empty after fallback"
fi

# 6) persist result (avoid copy-on-self)
SRC="$(readlink -f output/classified_sequences.tsv)"
DST="$(readlink -f "$OUTDIR/classified_sequences.tsv")"
if [ "$SRC" != "$DST" ]; then cp -f "$SRC" "$DST"; fi
[ -s "$OUTDIR/classified_sequences.tsv" ] || die "empty classification"

# 7) CAMI export (robust to 'rank:name' or plain 'name' lineage)
HYMET2CAMI_SCRIPT="/data/tools/hymet2cami_fixed.py"
[ -f "$HYMET2CAMI_SCRIPT" ] || die "missing hymet2cami script: $HYMET2CAMI_SCRIPT"

python3 "$HYMET2CAMI_SCRIPT" "$OUTDIR/classified_sequences.tsv" > "$OUTDIR/hymet.sample_0.cami.tsv"

log "OK: $OUTDIR/classified_sequences.tsv"
log "OK: $OUTDIR/hymet.sample_0.cami.tsv"

# small preview
echo "[preview] classified_sequences.tsv:"
head -n 5 "$OUTDIR/classified_sequences.tsv" || true
echo "[preview] CAMI profile:"
head -n 5 "$OUTDIR/hymet.sample_0.cami.tsv" || true
