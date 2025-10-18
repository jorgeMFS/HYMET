#!/usr/bin/env bash
# Run HYMET under progressive database ablations for the case-study dataset.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${CASE_ROOT}/manifest.tsv"
SAMPLE_ID=""
LEVELS="0,0.25,0.5,0.75,1.0"
TAXA=""
SEQMAP="${BENCH_ROOT}/db/ganon2/seqid2taxid.map"
BASE_FASTA="${HYMET_ROOT}/data/downloaded_genomes/combined_genomes.fasta"
OUT_ROOT="${CASE_ROOT}/ablation"
THREADS="${THREADS:-8}"
SEED=1337
MEASURE="${CASE_ROOT}/lib/measure.sh"
RUNTIME_TSV="${CASE_ROOT}/out/runtime_memory.tsv"

usage(){
  cat <<'USAGE'
Usage: run_ablation.sh --taxa TAXID1,TAXID2 [...] [--sample ID]
                       [--levels FRACTIONS] [--manifest TSV]
                       [--seqmap PATH] [--fasta PATH]
                       [--out DIR] [--threads N] [--seed N]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE_ID="$2"; shift 2;;
    --taxa) TAXA="$2"; shift 2;;
    --levels) LEVELS="$2"; shift 2;;
    --manifest) MANIFEST="$2"; shift 2;;
    --seqmap) SEQMAP="$2"; shift 2;;
    --fasta) BASE_FASTA="$2"; shift 2;;
    --out) OUT_ROOT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${TAXA}" ]] || die "Target TaxIDs must be provided via --taxa."

MANIFEST="$(resolve_path "${MANIFEST}")"
SEQMAP="$(resolve_path "${SEQMAP}")"
BASE_FASTA="$(resolve_path "${BASE_FASTA}")"
OUT_ROOT="$(resolve_path "${OUT_ROOT}")"
ensure_dir "${OUT_ROOT}"

[[ -s "${SEQMAP}" ]] || die "Sequence → taxid map missing: ${SEQMAP}"
[[ -s "${BASE_FASTA}" ]] || die "Reference FASTA missing: ${BASE_FASTA}"

# Determine sample (default first row in manifest)
if [[ -z "${SAMPLE_ID}" ]]; then
while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" || "${line}" == \#* || "${line}" == sample_id* ]] && continue
  IFS=$'\x1f' read -r SAMPLE_ID contigs _truthc _truthp _expected _citation <<<"$(manifest_split_line "${line}")"
  break
done < "${MANIFEST}"
  [[ -n "${SAMPLE_ID}" ]] || die "No sample rows found in manifest ${MANIFEST}"
fi

# Resolve contigs path for the selected sample
CONTIGS=""
TRUTH_CONTIGS=""
TRUTH_PROFILE=""
EXPECTED=""
CITATION=""
while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" || "${line}" == \#* || "${line}" == sample_id* ]] && continue
  IFS=$'\x1f' read -r sid contigs truth_contigs truth_profile expected citation <<<"$(manifest_split_line "${line}")"
  if [[ "${sid}" == "${SAMPLE_ID}" ]]; then
    CONTIGS="$(resolve_path "${contigs}")"
    TRUTH_CONTIGS="$(resolve_path "${truth_contigs}")"
    TRUTH_PROFILE="$(resolve_path "${truth_profile}")"
    break
  fi
done < "${MANIFEST}"

[[ -s "${CONTIGS}" ]] || die "Contigs for sample ${SAMPLE_ID} not found: ${CONTIGS}"

ABLATE_DIR="${OUT_ROOT}/refsets"
ensure_dir "${ABLATE_DIR}"

python3 "${CASE_ROOT}/ablate_db.py" \
  --fasta "${BASE_FASTA}" \
  --seqmap "${SEQMAP}" \
  --taxa "${TAXA}" \
  --levels "${LEVELS}" \
  --out-dir "${ABLATE_DIR}" \
  --prefix "combined_subset" \
  --seed "${SEED}"

SUMMARY_TSV="${OUT_ROOT}/ablation_summary.tsv"
if [[ ! -s "${SUMMARY_TSV}" ]]; then
  ensure_dir "$(dirname "${SUMMARY_TSV}")"
  cat <<'EOF' >"${SUMMARY_TSV}"
level_label	level_fraction	total_classified	assigned_species_pct	assigned_genus_pct	assigned_family_pct	assigned_higher_pct
EOF
fi

EVAL_SUMMARY="${OUT_ROOT}/ablation_eval_summary.tsv"
if [[ ! -s "${EVAL_SUMMARY}" ]]; then
  ensure_dir "$(dirname "${EVAL_SUMMARY}")"
  cat <<'EOF' >"${EVAL_SUMMARY}"
level_label	level_fraction	rank	F1	Precision	Recall	L1_total_variation_pctpts	BrayCurtis_pct	Contig_Accuracy_pct	Contig_Misassignment_pct
EOF
fi

backup="${BASE_FASTA}.case_backup"
if [[ ! -f "${backup}" ]]; then
  cp "${BASE_FASTA}" "${backup}"
fi

restore_fastas(){
  if [[ -f "${backup}" ]]; then
    mv -f "${backup}" "${BASE_FASTA}"
  fi
}
trap restore_fastas EXIT

for level_path in "${ABLATE_DIR}"/combined_subset.ablate*.fasta; do
  level_file="$(basename "${level_path}")"
  level_label="${level_file#combined_subset.ablate}"
  level_label="${level_label%.fasta}"
  level_fraction=$(python3 - <<'PY' "${level_label}"
import sys
label = sys.argv[1]
print(int(label)/100)
PY
)

  log "[ablation] Level ${level_label} (${level_fraction})"
  # Replace reference FASTA with the ablated version (except 000)
  if [[ "${level_label}" != "000" ]]; then
    cp "${level_path}" "${BASE_FASTA}"
  else
    cp "${backup}" "${BASE_FASTA}"
  fi

  level_out="${OUT_ROOT}/${SAMPLE_ID}/level_${level_label}"
  hymet_out="${level_out}/hymet"
  ensure_dir "${hymet_out}"

  THREADS="${THREADS}" "${MEASURE}" \
    --tool hymet \
    --sample "${SAMPLE_ID}_abl${level_label}" \
    --stage "ablation_${level_label}" \
    --out "${RUNTIME_TSV}" \
    -- "${BENCH_ROOT}/run_hymet.sh" \
         --sample "${SAMPLE_ID}" \
         --contigs "${CONTIGS}" \
         --out "${hymet_out}" \
         --threads "${THREADS}"

  profile="${hymet_out}/profile.cami.tsv"
  classified="${hymet_out}/classified_sequences.tsv"
  [[ -s "${classified}" ]] || { log "WARNING: classified_sequences.tsv missing for ${level_label}"; continue; }

  if [[ -n "${TRUTH_PROFILE}" && -s "${TRUTH_PROFILE}" && -n "${TRUTH_CONTIGS}" && -s "${TRUTH_CONTIGS}" ]]; then
    eval_stage="${SAMPLE_ID}_abl${level_label}"
    THREADS="${THREADS}" "${MEASURE}" \
      --tool hymet \
      --sample "${eval_stage}" \
      --stage "ablation_eval_${level_label}" \
      --out "${RUNTIME_TSV}" \
      -- "${BENCH_ROOT}/lib/run_eval.sh" \
            --sample "${SAMPLE_ID}" \
            --tool hymet \
            --pred-profile "${profile}" \
            --truth-profile "${TRUTH_PROFILE}" \
            --pred-contigs "${classified}" \
            --truth-contigs "${TRUTH_CONTIGS}" \
            --pred-fasta "${CONTIGS}" \
            --threads "${THREADS}" \
            --outdir "${hymet_out}/eval"

    profile_eval="${hymet_out}/eval/profile_summary.tsv"
    contig_eval="${hymet_out}/eval/contigs_per_rank.tsv"
    if [[ -s "${profile_eval}" ]]; then
      python3 - "${profile_eval}" "${contig_eval}" "${EVAL_SUMMARY}" "${level_label}" "${level_fraction}" <<'PY'
import csv, sys
from pathlib import Path
profile_path, contig_path, summary_path, level_label, level_fraction = sys.argv[1:6]
profile_rows = {}
with open(profile_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        profile_rows[row["rank"]] = row

contig_rows = {}
if contig_path and contig_path != "" and Path(contig_path).exists():
    with open(contig_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig_rows[row["rank"]] = row

with open(summary_path, "a", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    for rank, prow in profile_rows.items():
        contig_acc = contig_rows.get(rank, {}).get("accuracy_percent")
        if contig_acc is None or contig_acc == "":
            contig_acc = ""
            misassign = ""
        else:
            acc_val = float(contig_acc)
            contig_acc = f"{acc_val:.2f}"
            misassign = f"{max(0.0, 100.0 - acc_val):.2f}"
        writer.writerow([
            level_label,
            level_fraction,
            rank,
            prow.get("F1_%", ""),
            prow.get("Precision_%", ""),
            prow.get("Recall_%", ""),
            prow.get("L1_total_variation_pctpts", ""),
            prow.get("BrayCurtis_pct", ""),
            contig_acc,
            misassign,
        ])
PY
    fi
  fi

  python3 - "${SUMMARY_TSV}" "${classified}" "${level_label}" "${level_fraction}" <<'PY'
import csv, sys
summary_path, classified_path, level_label, level_fraction = sys.argv[1:5]
total = 0
by_rank = {"species": 0, "genus": 0, "family": 0, "order": 0, "class": 0, "phylum": 0, "superkingdom": 0}
with open(classified_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        total += 1
        rank = (row.get("Taxonomic Level") or "").strip().lower()
        if rank in by_rank:
            by_rank[rank] += 1
with open(summary_path, "a", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    if total == 0:
        writer.writerow([level_label, level_fraction, 0, 0, 0, 0, 0])
    else:
        species_pct = 100.0 * by_rank["species"] / total
        genus_pct = 100.0 * (by_rank["genus"] + by_rank["species"]) / total
        family_pct = 100.0 * (by_rank["family"] + by_rank["genus"] + by_rank["species"]) / total
        higher = 100.0 * (1.0 - (by_rank["species"] + by_rank["genus"] + by_rank["family"]) / total)
        writer.writerow([
            level_label,
            level_fraction,
            total,
            f"{species_pct:.2f}",
            f"{genus_pct:.2f}",
            f"{family_pct:.2f}",
            f"{higher:.2f}",
        ])
PY
done

restore_fastas
log "[ablation] Results written under ${OUT_ROOT}"

python3 "${CASE_ROOT}/plot_ablation.py" \
  --summary "${SUMMARY_TSV}" \
  --eval "${EVAL_SUMMARY}" \
  --outdir "${OUT_ROOT}/figures"
