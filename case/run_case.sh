#!/usr/bin/env bash
# Execute HYMET on real-world case-study samples defined in case/manifest.tsv.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${CASE_ROOT}/manifest.tsv"
OUT_ROOT="${CASE_ROOT}/out"
THREADS="${THREADS:-8}"
TOP_N="${TOP_N:-10}"
SANITY_METAPHLAN=0
METAPHLAN_CMD="${METAPHLAN_CMD:-metaphlan}"
METAPHLAN_OPTS="${METAPHLAN_OPTS:-}"
MEASURE="${CASE_ROOT}/lib/measure.sh"
RUNTIME_TSV="${OUT_ROOT}/runtime_memory.tsv"

usage(){
  cat <<'USAGE'
Usage: run_case.sh [--manifest TSV] [--out DIR] [--threads N]
                   [--top-n K] [--sanity-metaphlan]
                   [--metaphlan-cmd PATH] [--metaphlan-opts "..."]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --out) OUT_ROOT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --top-n) TOP_N="$2"; shift 2;;
    --sanity-metaphlan) SANITY_METAPHLAN=1; shift;;
    --metaphlan-cmd) METAPHLAN_CMD="$2"; shift 2;;
    --metaphlan-opts) METAPHLAN_OPTS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
OUT_ROOT="$(resolve_path "${OUT_ROOT}")"
ensure_dir "${OUT_ROOT}"

[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

append_summary_header(){
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    ensure_dir "$(dirname "${path}")"
    cat <<'EOF' >"${path}"
sample	metric	value
EOF
  fi
}

TOP_SUMMARY="${OUT_ROOT}/top_taxa_summary.tsv"
append_summary_header "${TOP_SUMMARY}"

while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" || "${line}" == \#* || "${line}" == sample_id* ]] && continue
  IFS=$'\0' read -r sample_id contigs expected_taxa citation <<<"$(manifest_split_line "${line}")"

  if [[ -z "${sample_id}" ]]; then
    log "Skipping manifest line with empty sample_id."
    continue
  fi

  contigs_abs="$(resolve_path "${contigs}")"
  if [[ ! -s "${contigs_abs}" ]]; then
    log "WARNING: sample ${sample_id} missing contigs at ${contigs_abs}; skipping."
    continue
  fi

  sample_out="${OUT_ROOT}/${sample_id}"
  hymet_out="${sample_out}/hymet"
  ensure_dir "${hymet_out}"

  log "[case] Running HYMET for ${sample_id}"
  THREADS="${THREADS}" "${MEASURE}" \
    --tool hymet \
    --sample "${sample_id}" \
    --stage run \
    --out "${RUNTIME_TSV}" \
    -- "${BENCH_ROOT}/run_hymet.sh" \
         --sample "${sample_id}" \
         --contigs "${contigs_abs}" \
         --out "${hymet_out}" \
         --threads "${THREADS}"

  profile="${hymet_out}/profile.cami.tsv"
  classified="${hymet_out}/classified_sequences.tsv"
  [[ -s "${profile}" ]] || die "HYMET profile missing for ${sample_id}: ${profile}"

  top_table="${sample_out}/top_taxa.tsv"
  python3 - "${profile}" "${top_table}" "${TOP_N}" <<'PY'
import csv, sys
profile_path, out_path, top_n = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = []
with open(profile_path) as fh:
    for line in fh:
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        taxid, rank, taxpath, taxpathsn, pct = line.split("\t")
        try:
            pct_val = float(pct)
        except ValueError:
            continue
        rows.append((pct_val, rank, taxid, taxpathsn, taxpath))
rows.sort(reverse=True)
with open(out_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Rank", "TaxID", "TaxPathSN", "TaxPath", "Percentage"])
    for pct, rank, taxid, taxpathsn, taxpath in rows[:top_n]:
        writer.writerow([rank, taxid, taxpathsn, taxpath, f"{pct:.6f}"])
PY

  python3 - "${TOP_SUMMARY}" "${sample_id}" "${top_table}" <<'PY'
import csv, sys, pathlib
summary_path, sample_id, table_path = sys.argv[1:4]
rows = []
with open(table_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        rows.append(row)
with open(summary_path, "a", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    for entry in rows:
        writer.writerow([sample_id, f"top_{entry['Rank']}", f"{entry['TaxPathSN']} ({entry['Percentage']})"])
PY

  cat > "${sample_out}/metadata.json" <<EOF
{
  "sample_id": "${sample_id}",
  "contigs": "${contigs_abs}",
  "expected_taxa": "${expected_taxa}",
  "citation": "${citation}",
  "hymet_profile": "${profile}",
  "hymet_classified": "${classified}"
}
EOF

  if [[ "${SANITY_METAPHLAN}" -eq 1 ]]; then
    metaphlan_out="${sample_out}/metaphlan"
    ensure_dir "${metaphlan_out}"
    mp_profile="${metaphlan_out}/profile.tsv"
    log "[case] Running MetaPhlAn sanity check for ${sample_id}"
    "${MEASURE}" \
      --tool metaphlan4 \
      --sample "${sample_id}" \
      --stage run \
      --out "${RUNTIME_TSV}" \
      -- "${METAPHLAN_CMD}" "${contigs_abs}" \
            --input_type fasta \
            --nproc "${THREADS}" \
            ${METAPHLAN_OPTS} \
            -o "${mp_profile}"

    if [[ -s "${mp_profile}" ]]; then
      comparison="${metaphlan_out}/comparison.tsv"
      python3 - "${profile}" "${mp_profile}" "${comparison}" <<'PY'
import csv, sys, math
hymet_path, mp_path, out_path = sys.argv[1:4]
def load_profile(path):
    prof = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            taxpathsn = parts[3]
            pct = float(parts[4])
            prof[taxpathsn] = pct
    return prof
hymet = load_profile(hymet_path)
mp = load_profile(mp_path)
taxa = sorted(set(hymet) | set(mp))
with open(out_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["TaxPathSN", "HYMET_Percent", "MetaPhlAn_Percent", "Absolute_Difference"])
    for taxon in taxa:
        h = hymet.get(taxon, 0.0)
        m = mp.get(taxon, 0.0)
        writer.writerow([taxon, f"{h:.6f}", f"{m:.6f}", f"{abs(h-m):.6f}"])
PY
    else
      log "WARNING: MetaPhlAn profile missing for ${sample_id}; comparison skipped."
    fi
  fi
done < "${MANIFEST}"

log "[case] Completed case-study run. Outputs in ${OUT_ROOT}"
