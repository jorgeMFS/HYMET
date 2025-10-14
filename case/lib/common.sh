#!/usr/bin/env bash
# Shared helpers for case-study scripts.

if [[ -n "${CASE_COMMON_SOURCED:-}" ]]; then
  return 0
fi
CASE_COMMON_SOURCED=1

set -o errexit
set -o nounset
set -o pipefail

_case__this="${BASH_SOURCE[0]}"
if [[ "${_case__this}" != */* ]]; then
  _case__this="./${_case__this}"
fi
_case__dir="$(cd "$(dirname "${_case__this}")" && pwd)"
export CASE_ROOT="$(cd "${_case__dir}/.." && pwd)"
export HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"
export BENCH_ROOT="$(cd "${CASE_ROOT}/../bench" && pwd)"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

ensure_dir(){
  local path="$1"
  mkdir -p "${path}"
}

resolve_path(){
  local input="$1"
  python3 - "$input" "$CASE_ROOT" <<'PY'
import os, sys
value, case_root = sys.argv[1], sys.argv[2]
if not value:
    print("", end="")
elif os.path.isabs(value):
    print(os.path.normpath(value), end="")
else:
    print(os.path.normpath(os.path.join(case_root, value)), end="")
PY
}

manifest_split_line(){
  local line="$1"
  python3 - <<'PY' -- "${line}"
import csv, sys
line = sys.argv[1]
row = next(csv.reader([line], delimiter='\t'))
sys.stdout.write("\0".join(row))
PY
}
