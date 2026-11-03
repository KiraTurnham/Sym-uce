#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=== [$(date)] Setting up Symbiodiniaceae LSU reference ==="

if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA is missing or empty: ${REF_FASTA}" >&2
  echo "Copy references/symbiodiniaceae_lsu.fasta to this path first." >&2
  exit 1
fi

mkdir -p "${REF_DIR}"
module purge
module load aligners/bowtie2/2.5.0

echo "Building Bowtie2 index at ${REF_INDEX_PREFIX}..."
bowtie2-build "${REF_FASTA}" "${REF_INDEX_PREFIX}"

echo "=== [$(date)] LSU reference index construction complete. ==="
