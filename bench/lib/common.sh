#!/usr/bin/env bash
# Shared helpers for bench scripts.

if [[ -n "${BENCH_COMMON_SOURCED:-}" ]]; then
  return 0
fi
BENCH_COMMON_SOURCED=1

set -o errexit
set -o nounset
set -o pipefail

_bench__this="${BASH_SOURCE[0]}"
if [[ "${_bench__this}" != */* ]]; then
  _bench__this="./${_bench__this}"
fi
_bench__dir="$(cd "$(dirname "${_bench__this}")" && pwd)"
export BENCH_ROOT="$(cd "${_bench__dir}/.." && pwd)"
export HYMET_ROOT="$(cd "${BENCH_ROOT}/.." && pwd)"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

ensure_dir(){
  local path="$1"
  mkdir -p "${path}"
}

resolve_path(){
  local input="$1"
  python3 - "$input" "$BENCH_ROOT" <<'PY'
import os, sys
value, bench_root = sys.argv[1], sys.argv[2]
if not value:
    print("", end="")
elif os.path.isabs(value):
    print(os.path.normpath(value), end="")
else:
    print(os.path.normpath(os.path.join(bench_root, value)), end="")
PY
}

need_cmd(){
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    log "ERROR: missing required command '${cmd}'"
    return 1
  fi
}

append_runtime_header(){
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    ensure_dir "$(dirname "${path}")"
    cat <<'EOF' >"${path}"
sample	tool	stage	wall_seconds	user_seconds	sys_seconds	max_rss_gb	io_input_mb	io_output_mb
EOF
  fi
}

manifest_split_line(){
  # Usage: manifest_split_line "<line>"
  local line="$1"
  python3 - <<'PY' -- "${line}"
import csv, sys
line = sys.argv[1]
row = next(csv.reader([line], delimiter="\t"))
sys.stdout.write("\0".join(row))
PY
}
