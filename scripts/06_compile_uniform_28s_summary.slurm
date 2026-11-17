#!/usr/bin/env bash
#SBATCH --job-name=sym_28s_summary
#SBATCH --output=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%j.out
#SBATCH --error=/scratch/%u/symbiodiniaceae_mining/slurm_logs/%x_%j.err
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/scripts"
if [[ ! -r "${SCRIPT_DIR}/config.sh" ]]; then
  echo "ERROR: Cannot find ${SCRIPT_DIR}/config.sh. Submit from the repository root." >&2
  exit 1
fi
source "${SCRIPT_DIR}/config.sh"

BLAST_DB="${REF_DIR}/sym_lsu_blastdb"
BLAST_OUT_DIR="${BASE_SCRATCH}/blast_results"
UNIFORM_OUT="${UNIFORM_OUT:-${BLAST_OUT_DIR}/all_175_colonies_uniform_read_and_contig_calls.tsv}"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Manifest is missing or empty: ${MANIFEST}" >&2
  exit 1
fi

DB_READY=false
for suffix in ndb nhr nin nal; do
  if [[ -s "${BLAST_DB}.${suffix}" ]]; then
    DB_READY=true
    break
  fi
done
if [[ "${DB_READY}" != true ]]; then
  echo "ERROR: BLAST database not found at ${BLAST_DB}; build it with makeblastdb first." >&2
  exit 1
fi

module purge
module load bio/blast/2.15.0+

mkdir -p "$(dirname "${UNIFORM_OUT}")" "${BLAST_OUT_DIR}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sym_28s_summary.XXXXXX")"
TMP_OUT="${UNIFORM_OUT}.tmp.$$"
trap 'rm -rf -- "${WORK_DIR}"; rm -f -- "${TMP_OUT}"' EXIT

printf "sample_id\tlsu_read_pairs\tr1_reads_with_hits\tdurusdinium_best_hits\tpct_durusdinium_of_hits\tnum_contigs\ttop_contig_hit\ttop_contig_pident\n" \
  > "${TMP_OUT}"

SAMPLE_COUNT=0
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

  R1="${BAM_DIR}/${SAMPLE}/${SAMPLE}_lsu_R1.fastq.gz"
  CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa"
  CONTIG_BLAST="${BLAST_OUT_DIR}/${SAMPLE}/${SAMPLE}_lsu_blastn.tsv"
  SAMPLE_WORK="${WORK_DIR}/${SAMPLE}"
  mkdir -p "${SAMPLE_WORK}"

  PAIRS=0
  READS_WITH_HITS=0
  DUR_HITS=0
  PCT_DUR="NA"
  if [[ -s "${R1}" ]]; then
    PAIRS=$(( $(gzip -cd "${R1}" | awk 'END { print NR }') / 4 ))
    if [[ "${PAIRS}" -gt 0 ]]; then
      gzip -cd "${R1}" \
        | awk 'NR%4==1{sub(/^@/, ">"); print} NR%4==2{print}' \
        | blastn -db "${BLAST_DB}" -task blastn -evalue 1e-10 \
            -max_target_seqs 5 -num_threads "${SLURM_CPUS_PER_TASK}" \
            -outfmt "6 qseqid stitle bitscore" \
        | sort -t $'\t' -k1,1 -k3,3gr \
        | awk -F'\t' '!seen[$1]++' > "${SAMPLE_WORK}/read_hits.tsv"

      READS_WITH_HITS=$(awk 'END { print NR+0 }' "${SAMPLE_WORK}/read_hits.tsv")
      if [[ "${READS_WITH_HITS}" -gt 0 ]]; then
        DUR_HITS=$(awk -F'\t' 'tolower($2) ~ /durusdinium/ { count++ } END { print count+0 }' \
          "${SAMPLE_WORK}/read_hits.tsv")
        PCT_DUR=$(awk -v d="${DUR_HITS}" -v b="${READS_WITH_HITS}" \
          'BEGIN { printf "%.2f", (d/b)*100 }')
      fi
    fi
  fi

  NUM_CONTIGS=0
  if [[ -s "${CONTIGS}" ]]; then
    NUM_CONTIGS=$(awk '/^>/ { count++ } END { print count+0 }' "${CONTIGS}")
  fi

  TOP_HIT="NA"
  TOP_PIDENT="NA"
  if [[ -s "${CONTIG_BLAST}" ]]; then
    TOP_CONTIG="$(sort -t $'\t' -k15,15gr "${CONTIG_BLAST}" | sed -n '1p')"
    if [[ -n "${TOP_CONTIG}" ]]; then
      TOP_HIT="$(printf '%s\n' "${TOP_CONTIG}" | cut -f3)"
      TOP_PIDENT="$(printf '%s\n' "${TOP_CONTIG}" | cut -f4)"
    fi
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" "${PAIRS}" "${READS_WITH_HITS}" "${DUR_HITS}" "${PCT_DUR}" \
    "${NUM_CONTIGS}" "${TOP_HIT}" "${TOP_PIDENT}" >> "${TMP_OUT}"
done < "${MANIFEST}"

if [[ "${SAMPLE_COUNT}" -lt 1 ]]; then
  echo "ERROR: No samples found in manifest: ${MANIFEST}" >&2
  exit 1
fi

mv -- "${TMP_OUT}" "${UNIFORM_OUT}"
echo "[$(date)] Wrote ${SAMPLE_COUNT} sample rows to ${UNIFORM_OUT}."

BACKUP_DIR="${HOME}/symbiodiniaceae_lsu_backup"
mkdir -p "${BACKUP_DIR}"
cp -- "${UNIFORM_OUT}" "${BACKUP_DIR}/"
if [[ -s "${REF_DIR}/symbiodiniaceae_lsu_id_map.tsv" ]]; then
  cp -- "${REF_DIR}/symbiodiniaceae_lsu_id_map.tsv" "${BACKUP_DIR}/"
fi
echo "[$(date)] Summary backed up to ${BACKUP_DIR}."
echo "Merged FASTQs were retained; remove them manually only after validating and archiving all required results."
