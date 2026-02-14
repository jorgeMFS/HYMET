#!/usr/bin/env bash
# =============================================================================
# Phase 5: Analyse HYMET results vs ZymoGut D6331 ground truth
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZYMOGUT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="$(cd "${ZYMOGUT_ROOT}/.." && pwd)"
HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"

WORK_DIR="${WORK_DIR:-${ZYMOGUT_ROOT}/work}"
FORCE=0

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) echo "Usage: 05_analyse.sh [--work-dir DIR] [--force]"; exit 0;;
    *) shift;;
  esac
done

TRUTH_DIR="${WORK_DIR}/truth"
HYMET_OUTPUT="${WORK_DIR}/hymet_output"
RESULTS_DIR="${WORK_DIR}/results"
FIGURES_DIR="${RESULTS_DIR}/figures"

TRUTH_CONTIGS="${TRUTH_DIR}/truth_contigs.tsv"
TRUTH_PROFILE="${TRUTH_DIR}/truth_profile.cami.tsv"
PREDICTED="${HYMET_OUTPUT}/classified_sequences.tsv"

log "ZymoGut D6331 Analysis"
log "======================"

[[ -s "${TRUTH_CONTIGS}" ]] || die "Truth contigs not found"
[[ -s "${PREDICTED}" ]] || die "HYMET output not found"

mkdir -p "${RESULTS_DIR}" "${FIGURES_DIR}"

log "Running analysis..."

# Step 1: Compute metrics and generate comparison tables
python3 - "${TRUTH_CONTIGS}" "${TRUTH_PROFILE}" "${PREDICTED}" "${RESULTS_DIR}" <<'PYEOF'
import sys, csv, json
from pathlib import Path
from collections import Counter
import numpy as np

truth_contigs_path = Path(sys.argv[1])
truth_profile_path = Path(sys.argv[2])
predicted_path = Path(sys.argv[3])
results_dir = Path(sys.argv[4])

print("Loading data...")

# Load truth contigs
truth_contigs = {}
with open(truth_contigs_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        truth_contigs[row['contig_id']] = {'species': row['species'], 'taxid': int(row['taxid']) if row['taxid'].isdigit() else 0}

# Load predictions
predictions = {}
with open(predicted_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cid = row['Query']
        lineage = row.get('Lineage', 'Unknown')
        ranks = {}
        if lineage and lineage != 'Unknown':
            for part in lineage.split(';'):
                if ':' in part:
                    rank, name = part.split(':', 1)
                    ranks[rank.strip().lower()] = name.strip()
        predictions[cid] = {'lineage': lineage, 'ranks': ranks}

# Load truth profile
truth_abundances = {}
with open(truth_profile_path) as f:
    for line in f:
        if line.startswith('#') or line.startswith('@'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 5 and parts[1] == 'species':
            species = parts[3].lower() if len(parts) > 3 else ''
            pct = float(parts[4]) if len(parts) > 4 else 0
            truth_abundances[species] = pct

print(f"  Truth contigs: {len(truth_contigs)}")
print(f"  Predictions: {len(predictions)}")

# Species comparison
predicted_species_counts = Counter()
for cid, pred in predictions.items():
    sp = pred['ranks'].get('species', 'Unclassified')
    predicted_species_counts[sp.lower()] += 1

total_contigs = len(predictions)
all_species = set(truth_abundances.keys()) | set(predicted_species_counts.keys())

# Comparison table
comparison_data = []
for sp in sorted(all_species):
    truth_pct = truth_abundances.get(sp, 0)
    pred_count = predicted_species_counts.get(sp, 0)
    pred_pct = (pred_count / total_contigs * 100) if total_contigs > 0 else 0
    comparison_data.append({'species': sp, 'truth_pct': truth_pct, 'predicted_count': pred_count, 'predicted_pct': pred_pct})

comparison_data.sort(key=lambda x: -x['truth_pct'])

with open(results_dir / 'comparison_table.tsv', 'w') as f:
    f.write("species\ttruth_pct\tpredicted_count\tpredicted_pct\tdiff_pct\n")
    for row in comparison_data:
        f.write(f"{row['species']}\t{row['truth_pct']:.4f}\t{row['predicted_count']}\t{row['predicted_pct']:.4f}\t{row['predicted_pct']-row['truth_pct']:.4f}\n")

# Profile metrics (species level)
species_list = sorted(all_species)
truth_vec = np.array([truth_abundances.get(sp, 0) / 100 for sp in species_list])
pred_vec = np.array([predicted_species_counts.get(sp, 0) / total_contigs if total_contigs > 0 else 0 for sp in species_list])

if truth_vec.sum() > 0:
    truth_vec = truth_vec / truth_vec.sum()
if pred_vec.sum() > 0:
    pred_vec = pred_vec / pred_vec.sum()

l1_distance = np.abs(truth_vec - pred_vec).sum()
bc_denom = (truth_vec + pred_vec).sum()
bray_curtis = np.abs(truth_vec - pred_vec).sum() / bc_denom if bc_denom > 0 else 0
correlation = np.corrcoef(truth_vec, pred_vec)[0, 1] if len(truth_vec) > 1 else 0

# Genus-level metrics (aggregate species to genus)
def to_genus(name):
    return name.split()[0] if name else name

truth_genus = {}
for sp, pct in truth_abundances.items():
    g = to_genus(sp)
    truth_genus[g] = truth_genus.get(g, 0) + pct

pred_genus_counts = Counter()
for cid, pred in predictions.items():
    g = pred['ranks'].get('genus', 'Unclassified')
    pred_genus_counts[g.lower()] += 1

all_genera = set(truth_genus.keys()) | set(pred_genus_counts.keys())
genus_list = sorted(all_genera)
truth_gvec = np.array([truth_genus.get(g, 0) / 100 for g in genus_list])
pred_gvec = np.array([pred_genus_counts.get(g, 0) / total_contigs if total_contigs > 0 else 0 for g in genus_list])

if truth_gvec.sum() > 0:
    truth_gvec = truth_gvec / truth_gvec.sum()
if pred_gvec.sum() > 0:
    pred_gvec = pred_gvec / pred_gvec.sum()

genus_l1 = np.abs(truth_gvec - pred_gvec).sum()
genus_bc_denom = (truth_gvec + pred_gvec).sum()
genus_bc = np.abs(truth_gvec - pred_gvec).sum() / genus_bc_denom if genus_bc_denom > 0 else 0
genus_corr = np.corrcoef(truth_gvec, pred_gvec)[0, 1] if len(truth_gvec) > 1 else 0

# Detection rates
truth_genera_set = {g for g, p in truth_genus.items() if p > 0}
pred_genera_set = {g for g, c in pred_genus_counts.items() if c > 0}
genus_recall = len(truth_genera_set & pred_genera_set) / len(truth_genera_set) if truth_genera_set else 0

# Genus comparison table
genus_comparison = []
for g in sorted(all_genera):
    t_pct = truth_genus.get(g, 0)
    p_count = pred_genus_counts.get(g, 0)
    p_pct = (p_count / total_contigs * 100) if total_contigs > 0 else 0
    genus_comparison.append({'genus': g, 'truth_pct': t_pct, 'predicted_count': p_count, 'predicted_pct': p_pct})

genus_comparison.sort(key=lambda x: -x['truth_pct'])

with open(results_dir / 'genus_comparison_table.tsv', 'w') as f:
    f.write("genus\ttruth_pct\tpredicted_count\tpredicted_pct\tdiff_pct\n")
    for row in genus_comparison:
        f.write(f"{row['genus']}\t{row['truth_pct']:.4f}\t{row['predicted_count']}\t{row['predicted_pct']:.4f}\t{row['predicted_pct']-row['truth_pct']:.4f}\n")

with open(results_dir / 'profile_metrics.tsv', 'w') as f:
    f.write("metric\tvalue\n")
    f.write(f"l1_distance\t{l1_distance:.6f}\n")
    f.write(f"bray_curtis\t{bray_curtis:.6f}\n")
    f.write(f"correlation\t{correlation:.6f}\n")
    f.write(f"genus_l1_distance\t{genus_l1:.6f}\n")
    f.write(f"genus_bray_curtis\t{genus_bc:.6f}\n")
    f.write(f"genus_correlation\t{genus_corr:.6f}\n")
    f.write(f"genus_recall\t{genus_recall:.6f}\n")
    f.write(f"total_contigs\t{total_contigs}\n")

print(f"\n  L1 distance: {l1_distance:.4f}")
print(f"  Bray-Curtis: {bray_curtis:.4f}")
print(f"  Correlation: {correlation:.4f}")
print(f"\n  Genus-level L1: {genus_l1:.4f}")
print(f"  Genus-level BC: {genus_bc:.4f}")
print(f"  Genus-level r:  {genus_corr:.4f}")
print(f"  Genus recall:   {genus_recall:.2%}")

# Summary JSON
with open(results_dir / 'analysis_summary.json', 'w') as f:
    json.dump({
        'sample': 'zymogut_d6331',
        'total_contigs': total_contigs,
        'metrics': {
            'l1_distance': float(l1_distance), 'bray_curtis': float(bray_curtis), 'correlation': float(correlation),
            'genus_l1_distance': float(genus_l1), 'genus_bray_curtis': float(genus_bc),
            'genus_correlation': float(genus_corr), 'genus_recall': float(genus_recall),
        }
    }, f, indent=2)

print("\nMetrics computation complete!")
PYEOF

log ""
log "Generating figures..."

# Step 2: Generate figures using dedicated plotting script (consistent with other case studies)
python3 "${ZYMOGUT_ROOT}/plot_zymogut.py" \
  --results-dir "${RESULTS_DIR}" \
  --figures-dir "${FIGURES_DIR}"

log ""
log "Analysis complete!"
log "Results in: ${RESULTS_DIR}"
log "Next step: Run 06_package.sh --publish"

cat > "${WORK_DIR}/status_05_analyse.json" <<EOF
{
  "phase": "05_analyse",
  "status": "complete",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "results_dir": "${RESULTS_DIR}"
}
EOF
