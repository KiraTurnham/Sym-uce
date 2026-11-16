#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 SAMPLE" >&2
  exit 2
fi

SAMPLE="$1"
if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
  echo "ERROR: Unsafe sample name: ${SAMPLE}" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ONE_SAMPLE_MANIFEST="$(mktemp "${TMPDIR:-/tmp}/sym_lsu_one_${SAMPLE}.XXXXXX")"
trap 'rm -f -- "${ONE_SAMPLE_MANIFEST}"' EXIT
printf '%s\n' "${SAMPLE}" > "${ONE_SAMPLE_MANIFEST}"

bash "${SCRIPT_DIR}/04_submit_lsu_new_samples.sh" "${ONE_SAMPLE_MANIFEST}"
