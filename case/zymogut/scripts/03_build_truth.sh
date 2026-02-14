#!/usr/bin/env bash
# =============================================================================
# Phase 3: Build ground truth files for ZymoGut D6331
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
THREADS="${THREADS:-8}"
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
Usage: 03_build_truth.sh [options]

Options:
  --work-dir DIR    Working directory
  --threads N       Number of threads
  --force           Force rebuild
  -h, --help        Show this message
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

ASSEMBLY_FASTA="${WORK_DIR}/assembly/assembly.fasta"
REFS_DIR="${WORK_DIR}/downloads/zymogut_refs"
TRUTH_DIR="${WORK_DIR}/truth"

log "ZymoGut D6331 Ground Truth Builder"
log "==================================="

[[ -s "${ASSEMBLY_FASTA}" ]] || die "Assembly not found: ${ASSEMBLY_FASTA}"
[[ -d "${REFS_DIR}" ]] || die "References not found: ${REFS_DIR}"

# Zymo D6331 zip extracts into genomes/ and ssrRNAs/ subdirectories.
# Use the genomes/ subdir when present, otherwise fall back to top-level.
if [[ -d "${REFS_DIR}/genomes" ]]; then
  GENOMES_DIR="${REFS_DIR}/genomes"
else
  GENOMES_DIR="${REFS_DIR}"
fi

mkdir -p "${TRUTH_DIR}"

# Step 1: Concatenate reference genomes
ALL_REFS="${TRUTH_DIR}/all_refs.fasta"
log "Step 1: Concatenating reference genomes..."
if [[ ! -s "${ALL_REFS}" || ${FORCE} -eq 1 ]]; then
  cat "${GENOMES_DIR}"/*.{fasta,fa,fna} 2>/dev/null > "${ALL_REFS}" || \
    cat "${GENOMES_DIR}"/*.fasta > "${ALL_REFS}" 2>/dev/null || \
    die "No FASTA files found in ${GENOMES_DIR}"
fi

# Step 2: Map contigs to references
PAF_FILE="${TRUTH_DIR}/contigs_vs_refs.paf"
log "Step 2: Mapping contigs to references..."
if [[ ! -s "${PAF_FILE}" || ${FORCE} -eq 1 ]]; then
  minimap2 -x asm5 -t "${THREADS}" --secondary=no "${ALL_REFS}" "${ASSEMBLY_FASTA}" > "${PAF_FILE}"
fi

# Step 3: Build truth files
TRUTH_CONTIGS="${TRUTH_DIR}/truth_contigs.tsv"
SEQID2TAXID="${TRUTH_DIR}/seqid2taxid.tsv"
TRUTH_PROFILE="${TRUTH_DIR}/truth_profile.cami.tsv"

log "Step 3: Building truth files..."

python3 - "${GENOMES_DIR}" "${PAF_FILE}" "${TRUTH_CONTIGS}" "${SEQID2TAXID}" "${TRUTH_PROFILE}" <<'PYEOF'
import sys, os, re, csv
from pathlib import Path
from collections import defaultdict, Counter

refs_dir = Path(sys.argv[1])
paf_file = Path(sys.argv[2])
truth_contigs_out = Path(sys.argv[3])
seqid2taxid_out = Path(sys.argv[4])
truth_profile_out = Path(sys.argv[5])

# ZymoGut D6331 species -> TaxID mapping
SPECIES_TO_TAXID = {
    "faecalibacterium prausnitzii": 853, "faecalibacterium": 853,
    "veillonella rogosae": 154447, "veillonella": 29465,
    "roseburia hominis": 301301, "roseburia": 841,
    "bacteroides fragilis": 817, "bacteroides": 816,
    "lactobacillus fermentum": 1613, "limosilactobacillus fermentum": 1613,
    "bifidobacterium adolescentis": 1680, "bifidobacterium": 1678,
    "fusobacterium nucleatum": 851, "fusobacterium": 848,
    "prevotella corporis": 28127, "prevotella": 838,
    "clostridioides difficile": 1496, "clostridium difficile": 1496,
    "akkermansia muciniphila": 239935, "akkermansia": 239934,
    "methanobrevibacter smithii": 2173, "methanobrevibacter": 2172,
    "salmonella enterica": 28901, "salmonella": 590,
    "escherichia coli": 562, "escherichia": 561,
    "enterococcus faecalis": 1351, "enterococcus": 1350,
    "saccharomyces cerevisiae": 4932, "saccharomyces": 4930,
    "candida albicans": 5476, "candida": 5475,
}

# Build ref_name -> taxid mapping
ref_to_species = {}
ref_to_taxid = {}

for fasta_file in refs_dir.glob("*"):
    if fasta_file.suffix.lower() not in ('.fasta', '.fa', '.fna'):
        continue
    filename = fasta_file.stem.lower().replace('_', ' ')
    filename_taxid = None
    for species, taxid in SPECIES_TO_TAXID.items():
        if species in filename or filename in species:
            filename_taxid = taxid
            filename_species = species
            break

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                ref_name = line[1:].split()[0]
                header = line.lower()
                header_taxid = None
                for species, taxid in sorted(SPECIES_TO_TAXID.items(), key=lambda x: -len(x[0])):
                    if species in header:
                        header_taxid = taxid
                        header_species = species
                        break
                if header_taxid:
                    ref_to_species[ref_name] = header_species
                    ref_to_taxid[ref_name] = header_taxid
                elif filename_taxid:
                    ref_to_species[ref_name] = filename_species
                    ref_to_taxid[ref_name] = filename_taxid
                else:
                    ref_to_species[ref_name] = filename
                    ref_to_taxid[ref_name] = 0

print(f"Mapped {len(ref_to_taxid)} reference sequences")

# Write seqid2taxid
with open(seqid2taxid_out, 'w') as f:
    f.write("seqid\ttaxid\tspecies\n")
    for ref_name in sorted(ref_to_taxid.keys()):
        f.write(f"{ref_name}\t{ref_to_taxid[ref_name]}\t{ref_to_species.get(ref_name, 'unknown')}\n")

# Parse PAF
contig_hits = defaultdict(list)
with open(paf_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue
        qname, qlen, qstart, qend = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        tname, mapq = parts[5], int(parts[11])
        coverage = (qend - qstart) / qlen if qlen > 0 else 0
        contig_hits[qname].append((tname, coverage, mapq))

# Best hit per contig
contig_truth = {}
for qname, hits in contig_hits.items():
    hits.sort(key=lambda x: (-x[1], -x[2]))
    best_ref, best_cov, best_mapq = hits[0]
    contig_truth[qname] = {
        'ref': best_ref,
        'species': ref_to_species.get(best_ref, 'unknown'),
        'taxid': ref_to_taxid.get(best_ref, 0),
        'coverage': best_cov,
        'mapq': best_mapq
    }

# Write contig truth
with open(truth_contigs_out, 'w') as f:
    f.write("contig_id\tspecies\ttaxid\tref_name\tcoverage\tmapq\n")
    for qname in sorted(contig_truth.keys()):
        t = contig_truth[qname]
        f.write(f"{qname}\t{t['species']}\t{t['taxid']}\t{t['ref']}\t{t['coverage']:.4f}\t{t['mapq']}\n")

print(f"Wrote contig truth: {len(contig_truth)} contigs")

# Write CAMI profile derived from assembly contig counts
species_counts = Counter()
for qname, t in contig_truth.items():
    species_counts[(t['species'], t['taxid'])] += 1

total_contigs = sum(species_counts.values())

with open(truth_profile_out, 'w') as f:
    f.write("#CAMI Submission for Taxonomic Profiling\n")
    f.write("@Version:0.9.1 @Ranks:superkingdom|phylum|class|order|family|genus|species @SampleID:zymogut_d6331_truth\n")
    f.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
    for (species, taxid), count in sorted(species_counts.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / total_contigs
        f.write(f"{taxid}\tspecies\t{taxid}\t{species}\t{pct:.6f}\n")

print(f"Wrote CAMI profile: {len(species_counts)} species (contig-derived)")
PYEOF

log "  Generated: ${TRUTH_CONTIGS}"
log "  Generated: ${SEQID2TAXID}"
log "  Generated: ${TRUTH_PROFILE}"

log ""
log "Ground truth building complete!"
log "Next step: Run 04_classify.sh"

CONTIG_COUNT=$(tail -n +2 "${TRUTH_CONTIGS}" | wc -l)
cat > "${WORK_DIR}/status_03_build_truth.json" <<EOF
{
  "phase": "03_build_truth",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "truth_contigs": "${TRUTH_CONTIGS}",
  "truth_profile": "${TRUTH_PROFILE}",
  "contig_count": ${CONTIG_COUNT}
}
EOF
