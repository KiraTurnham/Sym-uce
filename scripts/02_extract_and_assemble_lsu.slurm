#!/usr/bin/env bash
#SBATCH --job-name=sym_lsu_extract
#SBATCH --output=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%A_%a.out
#SBATCH --error=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%A_%a.err
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/scripts"
if [[ ! -r "${SCRIPT_DIR}/config.sh" ]]; then
  echo "ERROR: Cannot find ${SCRIPT_DIR}/config.sh." >&2
  echo "Submit this job from the repository root or set --chdir to it." >&2
  exit 1
fi
source "${SCRIPT_DIR}/config.sh"

if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: This script must run as a SLURM array task." >&2
  exit 1
fi
if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest file is missing or empty: ${MANIFEST}" >&2
  exit 1
fi
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA is missing or empty: ${REF_FASTA}" >&2
  exit 1
fi
if [[ ! -f "${REF_INDEX_PREFIX}.1.bt2" && ! -f "${REF_INDEX_PREFIX}.1.bt2l" ]]; then
  echo "ERROR: Bowtie2 index not found at ${REF_INDEX_PREFIX}" >&2
  echo "Run scripts/01_build_lsu_reference_index.sh first." >&2
  exit 1
fi
eval "$(mamba shell hook --shell bash)"
module purge
module load aligners/bowtie2/2.5.0
module load bio/samtools/1.19

mamba activate base

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample at SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}." >&2
  exit 1
fi
if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
  echo "ERROR: Unsafe sample name in manifest: ${SAMPLE}" >&2
  exit 1
fi

echo "Starting LSU read extraction for ${SAMPLE}"
echo "Node: $(hostname)"
echo "Date: $(date)"

R1="${FASTP_DIR}/${SAMPLE}${FASTQ_R1_SUFFIX}"
R2="${FASTP_DIR}/${SAMPLE}${FASTQ_R2_SUFFIX}"
if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "ERROR: Input FASTQ is missing or empty for ${SAMPLE}:" >&2
  echo "  R1: ${R1}" >&2
  echo "  R2: ${R2}" >&2
  exit 1
fi

SAMPLE_BAM="${BAM_DIR}/${SAMPLE}"
SAMPLE_ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE}"
mkdir -p "${SAMPLE_BAM}" "${ASSEMBLY_DIR}"

RAW_BAM="${SAMPLE_BAM}/${SAMPLE}_mapped_all.bam"
NAMESORT_BAM="${SAMPLE_BAM}/${SAMPLE}_namesorted.bam"
EXTRACTED_R1="${SAMPLE_BAM}/${SAMPLE}_lsu_R1.fastq.gz"
EXTRACTED_R2="${SAMPLE_BAM}/${SAMPLE}_lsu_R2.fastq.gz"
CONTIGS="${SAMPLE_ASSEMBLY}/final.contigs.fa"

for OUTPUT in \
  "${RAW_BAM}" "${NAMESORT_BAM}" "${EXTRACTED_R1}" "${EXTRACTED_R2}" \
  "${SAMPLE_ASSEMBLY}"; do
  if [[ -e "${OUTPUT}" ]]; then
    echo "ERROR: Output already exists for ${SAMPLE}: ${OUTPUT}" >&2
    echo "Move or remove the existing output before rerunning this sample." >&2
    exit 1
  fi
done

echo "[$(date)] Aligning ${SAMPLE} reads to LSU references..."
bowtie2 --very-sensitive-local -p "${THREADS}" \
  -x "${REF_INDEX_PREFIX}" -1 "${R1}" -2 "${R2}" \
  | samtools view -b -f 2 -o "${RAW_BAM}" -

echo "[$(date)] Sorting alignments by read name and extracting paired FASTQs..."
samtools sort -n -@ "${THREADS}" -o "${NAMESORT_BAM}" "${RAW_BAM}"
samtools fastq -@ "${THREADS}" \
  -1 "${EXTRACTED_R1}" -2 "${EXTRACTED_R2}" \
  -0 /dev/null -s /dev/null -n "${NAMESORT_BAM}"

rm -f -- "${RAW_BAM}" "${NAMESORT_BAM}"

MAPPED_READ_COUNT=$(( $(gzip -cd "${EXTRACTED_R1}" | wc -l) / 4 ))
echo "[$(date)] Recovered ${MAPPED_READ_COUNT} properly paired LSU-aligned read pairs."

if [[ "${MAPPED_READ_COUNT}" -lt "${MIN_MAPPED_PAIRS}" ]]; then
  echo "[$(date)] Fewer than ${MIN_MAPPED_PAIRS} LSU read pairs; skipping assembly."
  exit 0
fi

echo "[$(date)] Assembling target contigs with MEGAHIT..."
megahit -1 "${EXTRACTED_R1}" -2 "${EXTRACTED_R2}" \
  -t "${THREADS}" --min-contig-len "${MIN_CONTIG_LEN}" \
  -o "${SAMPLE_ASSEMBLY}"

if [[ -s "${CONTIGS}" ]]; then
  echo "[$(date)] Completed LSU extraction and assembly for ${SAMPLE}."
else
  echo "[$(date)] Warning: MEGAHIT produced no contigs of at least ${MIN_CONTIG_LEN} bp."
fi
