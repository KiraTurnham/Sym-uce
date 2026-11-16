#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 SAMPLE_MANIFEST" >&2
  echo "The manifest must contain only samples to submit, one sample name per line." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

INPUT_MANIFEST="$1"
if [[ ! -f "${INPUT_MANIFEST}" || ! -r "${INPUT_MANIFEST}" ]]; then
  echo "ERROR: Cannot read sample manifest: ${INPUT_MANIFEST}" >&2
  exit 1
fi
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA is missing or empty: ${REF_FASTA}" >&2
  echo "Copy references/symbiodiniaceae_lsu.fasta to this path first." >&2
  exit 1
fi

INDEX_READY=false
for INDEX_SUFFIX in bt2 bt2l; do
  if [[ -s "${REF_INDEX_PREFIX}.1.${INDEX_SUFFIX}" &&
        -s "${REF_INDEX_PREFIX}.2.${INDEX_SUFFIX}" &&
        -s "${REF_INDEX_PREFIX}.3.${INDEX_SUFFIX}" &&
        -s "${REF_INDEX_PREFIX}.4.${INDEX_SUFFIX}" &&
        -s "${REF_INDEX_PREFIX}.rev.1.${INDEX_SUFFIX}" &&
        -s "${REF_INDEX_PREFIX}.rev.2.${INDEX_SUFFIX}" ]]; then
    INDEX_READY=true
    break
  fi
done
if [[ "${INDEX_READY}" != true ]]; then
  echo "ERROR: Complete Bowtie2 index not found at ${REF_INDEX_PREFIX}." >&2
  echo "Run scripts/01_build_lsu_reference_index.sh first." >&2
  exit 1
fi
if [[ ! -d "${FASTP_DIR}" ]]; then
  echo "ERROR: fastp input directory does not exist: ${FASTP_DIR}" >&2
  exit 1
fi

mkdir -p "${BASE_SCRATCH}/sample_manifests" "${LOG_DIR}"
SAMPLE_COUNT=0
declare -A SEEN_SAMPLES=()
while IFS= read -r SAMPLE || [[ -n "${SAMPLE}" ]]; do
  if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: Sample manifest contains a blank line." >&2
    exit 1
  fi
  if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "ERROR: Unsafe sample name in manifest: ${SAMPLE}" >&2
    exit 1
  fi
  if [[ -n "${SEEN_SAMPLES[${SAMPLE}]:-}" ]]; then
    echo "ERROR: Duplicate sample in manifest: ${SAMPLE}" >&2
    exit 1
  fi
  SEEN_SAMPLES["${SAMPLE}"]=1

  R1="${FASTP_DIR}/${SAMPLE}${FASTQ_R1_SUFFIX}"
  R2="${FASTP_DIR}/${SAMPLE}${FASTQ_R2_SUFFIX}"
  if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "ERROR: Missing or empty fastp mate for ${SAMPLE}:" >&2
    echo "  R1: ${R1}" >&2
    echo "  R2: ${R2}" >&2
    exit 1
  fi
  SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done < "${INPUT_MANIFEST}"

if [[ "${SAMPLE_COUNT}" -lt 1 ]]; then
  echo "ERROR: Sample manifest is empty: ${INPUT_MANIFEST}" >&2
  exit 1
fi

SUBMIT_MANIFEST="${BASE_SCRATCH}/sample_manifests/submitted_$(date +%Y%m%d_%H%M%S)_$$.txt"
cp -- "${INPUT_MANIFEST}" "${SUBMIT_MANIFEST}"
cd "${REPO_DIR}"

echo "Submitting ${SAMPLE_COUNT} sample(s) from ${SUBMIT_MANIFEST}"
sbatch \
  --chdir="${REPO_DIR}" \
  --export="ALL,SAMPLE_MANIFEST=${SUBMIT_MANIFEST}" \
  --output="${LOG_DIR}/%x_%A_%a.out" \
  --error="${LOG_DIR}/%x_%A_%a.err" \
  --array="1-${SAMPLE_COUNT}%6" \
  "${SCRIPT_DIR}/02_extract_and_assemble_lsu.slurm"
