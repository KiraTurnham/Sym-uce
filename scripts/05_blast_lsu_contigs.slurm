#!/usr/bin/env bash
#SBATCH --job-name=sym_lsu_blast
#SBATCH --output=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%A_%a.out
#SBATCH --error=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%A_%a.err
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/scripts"
if [[ ! -r "${SCRIPT_DIR}/config.sh" ]]; then
  echo "ERROR: Cannot find ${SCRIPT_DIR}/config.sh. Submit from repo root." >&2
  exit 1
fi
source "${SCRIPT_DIR}/config.sh"

BLAST_DB="${REF_DIR}/sym_lsu_blastdb"
BLAST_OUT_DIR="${BASE_SCRATCH}/blast_results"

if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: Must run as a SLURM array task." >&2
  exit 1
fi
if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest missing or empty: ${MANIFEST}" >&2
  exit 1
fi
if [[ ! -f "${BLAST_DB}.nhr" && ! -f "${BLAST_DB}.nal" ]]; then
  echo "ERROR: BLAST database not found at ${BLAST_DB}. Run makeblastdb first." >&2
  exit 1
fi

module purge
module load bio/blast/2.15.0+

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
if [[ -z "${SAMPLE}" || ! "${SAMPLE}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
  echo "ERROR: Invalid or empty sample at line ${SLURM_ARRAY_TASK_ID}: ${SAMPLE}" >&2
  exit 1
fi

CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa"
SAMPLE_BLAST_DIR="${BLAST_OUT_DIR}/${SAMPLE}"
OUT_TSV="${SAMPLE_BLAST_DIR}/${SAMPLE}_lsu_blastn.tsv"

mkdir -p "${SAMPLE_BLAST_DIR}"

if [[ ! -s "${CONTIGS}" ]]; then
  echo "[$(date)] No assembled contigs found for ${SAMPLE} (${CONTIGS}); writing empty placeholder."
  : > "${OUT_TSV}"
  exit 0
fi

echo "[$(date)] Running blastn for ${SAMPLE} against ${BLAST_DB}..."

# Keep top 5 hits per contig to detect both primary Durusdinium trenchii and any co-occurring background genera
blastn \
  -query "${CONTIGS}" \
  -db "${BLAST_DB}" \
  -task blastn \
  -evalue 1e-20 \
  -max_target_seqs 5 \
  -num_threads "${SLURM_CPUS_PER_TASK}" \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qlen slen qstart qend sstart send evalue bitscore qcovs" \
  -out "${OUT_TSV}"

HIT_COUNT=$(wc -l < "${OUT_TSV}")
echo "[$(date)] Completed blastn for ${SAMPLE}: ${HIT_COUNT} alignment(s) written to ${OUT_TSV}."