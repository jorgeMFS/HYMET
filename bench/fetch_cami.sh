#!/usr/bin/env bash
# Idempotent helper to fetch CAMI assets referenced in the manifest.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${SCRIPT_DIR}/cami_manifest.tsv"
DRY_RUN=0
MAX_SAMPLES=0

usage(){
  cat <<'USAGE'
Usage: fetch_cami.sh [--manifest TSV] [--dry-run] [--max-samples N]

The manifest may optionally include columns with URLs named:
  contigs_url, truth_contigs_url, truth_profile_url
Only missing files with matching URLs are downloaded. Existing files are left untouched.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    --max-samples) MAX_SAMPLES="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
MANIFEST_DIR="$(dirname "${MANIFEST}")"
[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

CACHE_DIR="${MANIFEST_DIR}/tmp_downloads"
ensure_dir "${CACHE_DIR}"

ARCHIVE_FETCHED=()
DERIVED_GENERATED=0

SAMPLE0_ARCHIVE="https://frl.publisso.de/data/frl:6421672/dataset/2017.12.29_11.37.26_sample_0_contigs.tar"
SAMPLE0_CONTIG_MEMBER="2017.12.29_11.37.26_sample_0/contigs/anonymous_gsa.fasta.gz"
SAMPLE0_MAPPING_MEMBER="2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv.gz"
SAMPLE0_PROFILE_MEMBER="taxonomic_profile_0.txt"
SAMPLE0_BASE_DIR="/data/cami/sample_0"
SAMPLE0_FNA="/data/cami/sample_0.fna"
SAMPLE0_MAPPING="${SAMPLE0_BASE_DIR}/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv"

cache_archive_path(){
  local url="$1"
  local hash
  hash="$(printf '%s' "${url}" | sha1sum | cut -d' ' -f1)"
  printf '%s/%s.tar' "${CACHE_DIR}" "${hash}"
}

download_to_file(){
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    log "using cached ${dest}"
    return 0
  fi
  ensure_dir "$(dirname "${dest}")"
  log "downloading ${url}"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 4 -s 4 -o "${dest}" "${url}" || die "aria2c failed for ${url}"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "${dest}" "${url}" || die "curl failed for ${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${dest}" "${url}" || die "wget failed for ${url}"
  else
    die "No downloader (aria2c/curl/wget) available for ${url}"
  fi
}

extract_from_archive(){
  local archive="$1"
  local member="$2"
  local dest="$3"
  ensure_dir "$(dirname "${dest}")"
  if [[ "${member}" == *.gz && "${dest}" != *.gz ]]; then
    tar -xOf "${archive}" "${member}" | gzip -dc > "${dest}" || die "Failed to extract ${member} from ${archive}"
  else
    tar -xOf "${archive}" "${member}" > "${dest}" || die "Failed to extract ${member} from ${archive}"
  fi
}

ensure_sample0_assets(){
  if [[ -s "${SAMPLE0_FNA}" && -s "${SAMPLE0_MAPPING}" ]]; then
    return 0
  fi
  if [[ ${DRY_RUN} -eq 1 ]]; then
    log "would extract sample 0 assets into /data/cami"
    return 0
  fi
  local archive_path
  archive_path="$(cache_archive_path "${SAMPLE0_ARCHIVE}")"
  download_to_file "${SAMPLE0_ARCHIVE}" "${archive_path}"
  log "extracting sample 0 contigs to ${SAMPLE0_FNA}"
  extract_from_archive "${archive_path}" "${SAMPLE0_CONTIG_MEMBER}" "${SAMPLE0_FNA}"
  log "extracting sample 0 mapping to ${SAMPLE0_MAPPING}"
  extract_from_archive "${archive_path}" "${SAMPLE0_MAPPING_MEMBER}" "${SAMPLE0_MAPPING}"
  local profile_dest="/data/cami/sample_0/taxonomic_profile_0.txt"
  extract_from_archive "${archive_path}" "${SAMPLE0_PROFILE_MEMBER}" "${profile_dest}"
}

generate_derived_samples(){
  if [[ ${DERIVED_GENERATED} -eq 1 ]]; then
    return 0
  fi
  ensure_sample0_assets
  if [[ ${DRY_RUN} -eq 1 ]]; then
    log "would generate derived CAMI subsets under ${MANIFEST_DIR}/data"
    DERIVED_GENERATED=1
    return 0
  fi
  local script_path="${SCRIPT_DIR%/*}/tools/generate_cami_subsets.py"
  python3 "${script_path}" \
    --fasta "${SAMPLE0_FNA}" \
    --mapping "${SAMPLE0_MAPPING}" \
    --outdir "${MANIFEST_DIR}/data" || die "generate_cami_subsets.py failed"
  DERIVED_GENERATED=1
}

handle_missing_without_url(){
  local abs_path="$1"
  local rel_path="$2"
  local sample_name
  sample_name="$(basename "$(dirname "${rel_path}")")"
  case "${sample_name}" in
    cami_i_lc|cami_i_mc|cami_i_hc|cami_ii_mousegut|cami_ii_marine|cami_ii_strainmadness)
      generate_derived_samples
      if [[ -s "${abs_path}" ]]; then
        log "generated: ${abs_path}"
        return 0
      fi
      return 1
      ;;
    *)
      return 1
      ;;
  esac
}

python3 - "$MANIFEST" "$MAX_SAMPLES" <<'PY' | while IFS=$'	' read -r path url; do
import csv, sys
manifest = sys.argv[1]
max_samples = int(sys.argv[2]) if sys.argv[2] and sys.argv[2] != "0" else None
rows = []
with open(manifest, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for idx, row in enumerate(reader):
        if max_samples is not None and idx >= max_samples:
            break
        paths = {
            "contigs_fa": row.get("contigs_fa", ""),
            "truth_contigs_tsv": row.get("truth_contigs_tsv", ""),
            "truth_profile_tsv": row.get("truth_profile_tsv", ""),
        }
        urls = {
            "contigs_fa": row.get("contigs_url", ""),
            "truth_contigs_tsv": row.get("truth_contigs_url", ""),
            "truth_profile_tsv": row.get("truth_profile_url", ""),
        }
        for key, path in paths.items():
            if not path:
                continue
            url = urls.get(key) or ""
            print(f"{path}\t{url}")
PY
  [[ -z "${path}" ]] && continue
  abs_path="$(resolve_path "${path}" "${MANIFEST_DIR}")"
  url="${url}"
  if [[ -s "${abs_path}" ]]; then
    log "exists: ${abs_path}"
    continue
  fi
  if [[ -z "${url}" ]]; then
    rel_path="${path}"
    if handle_missing_without_url "${abs_path}" "${rel_path}"; then
      continue
    fi
    log "missing and no URL: ${abs_path}"
    continue
  fi
  archive_member=""
  archive_url="${url}"
  if [[ "${url}" == *"#"* ]]; then
    archive_member="${url#*#}"
    archive_url="${url%%#*}"
  fi
  if [[ -n "${archive_member}" ]]; then
    if [[ ${DRY_RUN} -eq 1 ]]; then
      log "would extract ${archive_url}#${archive_member} -> ${abs_path}"
      continue
    fi
    archive_path="$(cache_archive_path "${archive_url}")"
    download_to_file "${archive_url}" "${archive_path}"
    log "extracting ${archive_member} -> ${abs_path}"
    extract_from_archive "${archive_path}" "${archive_member}" "${abs_path}"
    continue
  fi
  ensure_dir "$(dirname "${abs_path}")"
  if [[ ${DRY_RUN} -eq 1 ]]; then
    log "would download ${archive_url} -> ${abs_path}"
    continue
  fi
  download_to_file "${archive_url}" "${abs_path}"
done
