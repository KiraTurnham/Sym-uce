#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 SAMPLE_MANIFEST" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

MANIFEST="$1"
OUT_DIR="${BASE_SCRATCH}/geneious_export"
CONTIGS_OUT="${OUT_DIR}/all_161_colonies_28S_LSU_contigs.fasta"
TREE_OUT="${OUT_DIR}/all_161_colonies_plus_LaJeunesse2018_refs.fasta"
CONTIGS_TMP="${CONTIGS_OUT}.tmp.$$"
TREE_TMP="${TREE_OUT}.tmp.$$"
BACKUP_DIR="${HOME}/symbiodiniaceae_lsu_backup/geneious_export"
BACKUP_CONTIGS="${BACKUP_DIR}/$(basename "${CONTIGS_OUT}")"
BACKUP_TREE="${BACKUP_DIR}/$(basename "${TREE_OUT}")"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest is missing or empty: ${MANIFEST}" >&2
  exit 1
fi
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA is missing or empty: ${REF_FASTA}" >&2
  exit 1
fi
for OUTPUT in "${CONTIGS_OUT}" "${TREE_OUT}" "${BACKUP_CONTIGS}" "${BACKUP_TREE}"; do
  if [[ -e "${OUTPUT}" ]]; then
    echo "ERROR: Output already exists: ${OUTPUT}" >&2
    echo "Move or remove the specific existing file before rerunning." >&2
    exit 1
  fi
done

mkdir -p "${OUT_DIR}" "${BACKUP_DIR}"
trap 'rm -f -- "${CONTIGS_TMP}" "${TREE_TMP}"' EXIT
: > "${CONTIGS_TMP}"

SAMPLE_COUNT=0
CONTIG_COUNT=0
declare -A SEEN_SAMPLES=()
while IFS= read -r SAMPLE || [[ -n "${SAMPLE}" ]]; do
  if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "ERROR: Invalid or empty sample ID in manifest: ${SAMPLE}" >&2
    exit 1
  fi
  if [[ -n "${SEEN_SAMPLES[${SAMPLE}]:-}" ]]; then
    echo "ERROR: Duplicate sample in manifest: ${SAMPLE}" >&2
    exit 1
  fi
  SEEN_SAMPLES["${SAMPLE}"]=1
  SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

  CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa"
  if [[ ! -s "${CONTIGS}" ]]; then
    echo "No contigs for ${SAMPLE}; skipping from Geneious FASTA."
    continue
  fi

  SAMPLE_CONTIGS=$(awk -v sample="${SAMPLE}" '
    /^>/ {
      n++
      len = "NA"
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^len=/) {
          sub(/^len=/, "", $i)
          len = $i
        }
      }
      print ">" sample "|contig_" n "|len_" len
      next
    }
    { print }
  ' "${CONTIGS}")
  if [[ -z "${SAMPLE_CONTIGS}" ]]; then
    echo "ERROR: No FASTA records found in ${CONTIGS}" >&2
    exit 1
  fi
  printf '%s\n' "${SAMPLE_CONTIGS}" >> "${CONTIGS_TMP}"
  SAMPLE_CONTIG_COUNT=$(awk '/^>/ { count++ } END { print count+0 }' "${CONTIGS}")
  CONTIG_COUNT=$((CONTIG_COUNT + SAMPLE_CONTIG_COUNT))
done < "${MANIFEST}"

if [[ "${SAMPLE_COUNT}" -lt 1 || "${CONTIG_COUNT}" -lt 1 ]]; then
  echo "ERROR: Manifest had no samples or no colony contigs to export." >&2
  exit 1
fi

cat "${CONTIGS_TMP}" "${REF_FASTA}" > "${TREE_TMP}"
mv -- "${CONTIGS_TMP}" "${CONTIGS_OUT}"
mv -- "${TREE_TMP}" "${TREE_OUT}"
cp -- "${CONTIGS_OUT}" "${TREE_OUT}" "${BACKUP_DIR}/"

TREE_COUNT=$(awk '/^>/ { count++ } END { print count+0 }' "${TREE_OUT}")
echo "Exported ${CONTIG_COUNT} colony contigs to ${CONTIGS_OUT}."
echo "Exported ${TREE_COUNT} total sequences to ${TREE_OUT}."
echo "Backed up FASTA exports to ${BACKUP_DIR}."
