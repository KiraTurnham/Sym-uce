# Symbiodiniaceae LSU Read Extraction

This Sedna SLURM workflow recruits paired-end reads to a Symbiodiniaceae
large-subunit ribosomal RNA (LSU/28S) reference, extracts properly paired
reads, assembles them with MEGAHIT, and can classify contigs with BLASTn.
Lane-level fastp outputs from the same colony must be concatenated by mate
before the extraction/assembly step so that each manifest entry represents one
colony, not one sequencing lane.

The LSU panel is an enrichment and comparison reference, not a comprehensive
taxonomic database. Neither Bowtie2 recruitment nor a top BLAST hit alone
establishes a species identification. Closely related templates may share
similar LSU sequences, and samples may contain multiple symbionts.

## Requirements on Sedna

- NOAA Fisheries Sedna SLURM and Environment Modules. The recorded cluster
  version is SLURM `25.11.6`; SLURM is site-managed, not installed by this
  repository.
- `bio/fastp/0.23.2` for the upstream trimming step (not run by these scripts).
- `aligners/bowtie2/2.5.0` and `bio/samtools/1.19` for recruitment/extraction.
- Mamba environment `mg-assembly` containing MEGAHIT.
- `bio/blast/2.15.0+` for BLAST database construction and BLASTn.
- `R/4.4.1` for downstream R analyses; no R script is included here.
- Bash, Python 3 (standard library), `gzip`, and standard GNU utilities.

Module and environment versions below are the versions reported for the Sedna
analysis. Verify them on the cluster and record executable versions with the
run logs; module names alone do not guarantee package patch versions.

```bash
module purge
module load bio/fastp/0.23.2
fastp --version

module purge
module load aligners/bowtie2/2.5.0
module load bio/samtools/1.19
bowtie2 --version
samtools --version

export MEGAHIT_ENV="${MEGAHIT_ENV:-mg-assembly}"
eval "$(mamba shell hook --shell bash)"
mamba activate "${MEGAHIT_ENV}"
megahit --version

module purge
module load bio/blast/2.15.0+
blastn -version
makeblastdb -version

module purge
module load R/4.4.1
R --version
sbatch --version
```

The extraction worker runs non-interactively and therefore initializes Mamba
with `eval "$(mamba shell hook --shell bash)"` before activating
`MEGAHIT_ENV`. The same hook is required in any other `sbatch` script that
activates a Mamba environment. The extraction worker loads Bowtie2 and SAMtools
modules itself; the BLAST worker loads BLAST+ itself.
`MEGAHIT_ENV` defaults to `mg-assembly` in `scripts/config.sh` and can be
overridden in the submitting shell.

## Setup

Clone the repository and create the scratch directories. Keep I/O-intensive
inputs, intermediate files, and outputs on `/scratch`, not NFS-backed home
storage.

```bash
git clone https://github.com/KiraTurnham/Sym-uce.git
cd Sym-uce

export SYMBIO_BASE_DIR="/scratch/${USER}/symbiodiniaceae_mining"
export FASTP_INPUT_DIR="${SYMBIO_BASE_DIR}/merged_fastp"
mkdir -p "${SYMBIO_BASE_DIR}"/{ref,sample_manifests,mapped_bam,assemblies,slurm_logs,blast_results,merged_fastp}
cp references/symbiodiniaceae_lsu.fasta \
  "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu.fasta"
```

The expected Sedna layout is:

```text
/scratch/kturnham/
├── Sym-uce/                                      # cloned repository; run scripts here
│   ├── README_sym.md
│   ├── references/
│   └── scripts/
├── symbiodiniaceae_mining/                       # SYMBIO_BASE_DIR: indexes, manifests, results
│   ├── ref/
│   │   ├── symbiodiniaceae_lsu.fasta
│   │   ├── symbiodiniaceae_lsu_blast_clean.fasta
│   │   ├── symbiodiniaceae_lsu_id_map.tsv
│   │   ├── sym_lsu_bt2.*                         # Bowtie2 index files
│   │   └── sym_lsu_blastdb.*                     # BLAST database files
│   ├── sample_manifests/
│   ├── mapped_bam/
│   ├── assemblies/
│   ├── blast_results/
│   ├── merged_fastp/                             # temporary colony-level merged FASTQs
│   └── slurm_logs/
├── phyluce_project/
│   ├── workflow_fast_megahit/
│   │   └── fastp_trimmed/                        # original lane-level fastp output; do not modify
│   │       ├── <sample>_<lane>_R1_trimmed.fastq.gz
│   │       └── <sample>_<lane>_R2_trimmed.fastq.gz
│   └── qc_reports/
│       └── trimmed_multiqc_new_only/
│           └── new_samples_only_manifest.txt    # derived from four project demultiplex tables
```

The checked-in scripts run from `/scratch/${USER}/Sym-uce`. They read the
colony-level merged FASTQs from `merged_fastp/` and write indexes, manifests,
alignments, assemblies, and logs under `symbiodiniaceae_mining/`.

The reference contains 75 LSU sequences from the Symbiodiniaceae family in
the LaJeunesse et al. systematic revision alignment, plus two 28S records for
*Durusdinium*. Its seven represented genera include *Cladocopium*, so reads
from an unrepresented species such as a potential *Cladocopium* C40 can still
be recruited by genus-level LSU similarity. The reference does not establish
species identity. Provenance and limitations are documented in
[`references/README.md`](references/README.md) and
[`references/lsu_reference_sources.tsv`](references/lsu_reference_sources.tsv).

Raw, lane-level outputs are expected to follow names such as
`1_icra5_S18_L005_R1_trimmed.fastq.gz` and the matching R2 filename. They were
trimmed upstream with fastp `0.23.2`; this repository does not rerun fastp.
Keep the original `fastp_trimmed/` directory intact.

Build the Bowtie2 index on Sedna:

```bash
bash scripts/01_build_lsu_reference_index.sh
```

## Sample manifest provenance and lane merging

Use the *new-run-only* manifest, not a directory scan: older sequencing runs
also occur in `fastp_trimmed/`. The reported `new_samples_only_manifest.txt`
contains 175 colony IDs and was derived from the four NovaSeq X Plus
Demultiplex_Stats.csv lane tables for Genohub Project `#4620737`, excluding
legacy runs present in the shared fastp directory. The four source tables,
their checksums, and the exact extraction/validation command are not included
in this repository; retain them with the analysis archive to make that
filtering independently reproducible. The manifest must contain one colony
ID per line, without lane or read suffixes.

Copy the verified manifest to scratch, then concatenate the matching lane
FASTQs per colony. Concatenating gzip streams preserves valid gzip input and
does not decompress/recompress the reads. Run the I/O-heavy step on a compute
node, not the login node:

```bash
export SRC="/scratch/${USER}/phyluce_project/workflow_fast_megahit/fastp_trimmed"
export MERGED="${SYMBIO_BASE_DIR}/merged_fastp"
export MANIFEST="/scratch/${USER}/phyluce_project/qc_reports/trimmed_multiqc_new_only/new_samples_only_manifest.txt"
cp "${MANIFEST}" "${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt"
export MANIFEST="${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt"

srun -p standard -c 1 --mem=2G --time=00:15:00 bash -c '
  set -euo pipefail
  shopt -s nullglob
  sample_count=0
  trap '\''rm -f -- "${tmp1:-}" "${tmp2:-}"'\'' EXIT
  mkdir -p "${MERGED}"
  while IFS= read -r sample || [[ -n "${sample}" ]]; do
    [[ "${sample}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
      echo "Invalid sample ID: ${sample}" >&2; exit 1;
    }
    r1=( "${SRC}/${sample}_"*_R1_trimmed.fastq.gz )
    r2=( "${SRC}/${sample}_"*_R2_trimmed.fastq.gz )
    [[ "${#r1[@]}" -gt 0 && "${#r1[@]}" -eq "${#r2[@]}" ]] || {
      echo "Missing/mismatched lane mates for ${sample}" >&2; exit 1;
    }
    out1="${MERGED}/${sample}_R1_trimmed.fastq.gz"
    out2="${MERGED}/${sample}_R2_trimmed.fastq.gz"
    [[ ! -e "${out1}" && ! -e "${out2}" ]] || {
      echo "Merged output already exists for ${sample}" >&2; exit 1;
    }
    tmp1="${out1}.tmp.$$"
    tmp2="${out2}.tmp.$$"
    for i in "${!r1[@]}"; do
      expected_r2="${r1[$i]/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
      [[ "${r2[$i]}" == "${expected_r2}" ]] || {
        echo "R1/R2 lane mismatch for ${sample}: ${r1[$i]}" >&2; exit 1;
      }
    done
    cat "${r1[@]}" > "${tmp1}"
    cat "${r2[@]}" > "${tmp2}"
    gzip -t "${tmp1}" "${tmp2}"
    mv -- "${tmp1}" "${out1}"
    mv -- "${tmp2}" "${out2}"
    sample_count=$((sample_count + 1))
  done < "${MANIFEST}"
  [[ "${sample_count}" -gt 0 ]] || {
    echo "Sample manifest is empty" >&2; exit 1;
  }
'
```

Set `FASTP_INPUT_DIR` to `merged_fastp/` in the same shell before submitting
the extraction jobs. The merge step fails rather than overwriting existing
merged FASTQs; inspect/remove only the specific stale merged files before
restarting a failed merge. Preserve the source fastp files.

## LSU reference index, extraction, and assembly

Build the Bowtie2 index once after copying the reference:

```bash
bash scripts/01_build_lsu_reference_index.sh
```

The workflow uses 8 CPUs, 4 GB RAM, and 30 minutes per sample on `-p standard`.
The array is throttled to six simultaneous tasks (`--array=1-N%6`). This
resource request is based on the reported `seff` measurements; monitor the
full batch and adjust if any task approaches its time or memory limit. The
submitter validates every sample and both FASTQ mates, snapshots the manifest,
creates the log directory, and submits the array:

```bash
export FASTP_INPUT_DIR="${SYMBIO_BASE_DIR}/merged_fastp"
bash scripts/03_submit_lsu_one.sh SAMPLE_ID
# Review the pilot, then exclude its ID from the full verified cohort manifest:
grep -vxF 'PILOT_SAMPLE_ID' \
  "${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt" \
  > "${SYMBIO_BASE_DIR}/sample_manifests/remaining_new_samples.txt"
# Submit only the remaining colony IDs:
bash scripts/04_submit_lsu_new_samples.sh \
  "${SYMBIO_BASE_DIR}/sample_manifests/remaining_new_samples.txt"
```

The remaining manifest must be created from the verified 175-colony manifest
and have the successful pilot ID removed. Do not derive it by scanning FASTQ
filenames or include previous runs. The worker extracts paired reads for every
sample and assembles only samples with at least `MIN_MAPPED_PAIRS=10`; MEGAHIT
uses `MIN_CONTIG_LEN=300`. It refuses to overwrite existing per-sample
outputs. Retain the merged FASTQs until all extraction/assembly jobs have
finished and been checked; then remove only the merged FASTQs to conserve
scratch quota (the original lane-level fastp files remain untouched):

```bash
rm -- "${SYMBIO_BASE_DIR}"/merged_fastp/*_R[12]_trimmed.fastq.gz
```

Outputs include `mapped_bam/<SAMPLE>/<SAMPLE>_lsu_R1.fastq.gz` and R2,
`assemblies/<SAMPLE>/final.contigs.fa` when MEGAHIT runs, and SLURM logs under
`slurm_logs/`. The extracted R1/R2 files remain available for downstream
read-level analyses, including low-coverage samples skipped for assembly.

## BLAST reference and contig classification

The reference has 77 sequences: 75 Symbiodiniaceae-family records from the
LaJeunesse et al. (2018) LSU alignment and two additional *Durusdinium* 28S
records (`KF740689.1`, `KY131785.1`). `makeblastdb -parse_seqids` requires
short identifiers; long local IDs can exceed BLAST's 50-character identifier
limit. Prepare deterministic short `lcl|symLSU####` IDs while preserving the
original FASTA deflines in the description and in a TSV ID map:

```bash
python3 scripts/prepare_lsu_blast_db.py \
  --input "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu.fasta" \
  --fasta-out "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu_blast_clean.fasta" \
  --id-map-out "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu_id_map.tsv"

module purge
module load bio/blast/2.15.0+
makeblastdb \
  -in "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu_blast_clean.fasta" \
  -dbtype nucl -parse_seqids \
  -title "Symbiodiniaceae_LSU_LaJeunesse2018_Durusdinium" \
  -out "${SYMBIO_BASE_DIR}/ref/sym_lsu_blastdb"
```

The original FASTA header can be recovered by looking up the `sseqid` in
`symbiodiniaceae_lsu_id_map.tsv`; BLAST output also includes `stitle`. Rebuild
the database from the cleaned FASTA if an earlier database was made directly
from the original long-header FASTA. The database and generated clean FASTA/map
are scratch outputs and should be archived with the analysis metadata.

[`scripts/05_blast_lsu_contigs.slurm`](scripts/05_blast_lsu_contigs.slurm)
queries each `assemblies/<SAMPLE>/final.contigs.fa` against this database with
BLAST+ `2.15.0+`, `blastn -evalue 1e-20 -max_target_seqs 5`. It writes one TSV
per sample under `blast_results/<SAMPLE>/`. This script is an array worker;
submit it with the same verified cohort manifest and an array range matching
its line count:

```bash
BLAST_MANIFEST="${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt"
[[ -s "${BLAST_MANIFEST}" ]] || { echo "Missing/empty BLAST manifest" >&2; exit 1; }
BLAST_N=$(awk 'END { print NR }' "${BLAST_MANIFEST}")
[[ "${BLAST_N}" -gt 0 ]] || { echo "BLAST manifest is empty" >&2; exit 1; }
sbatch --chdir="${PWD}" \
  --export="ALL,SAMPLE_MANIFEST=${BLAST_MANIFEST}" \
  --output="${SYMBIO_BASE_DIR}/slurm_logs/%x_%A_%a.out" \
  --error="${SYMBIO_BASE_DIR}/slurm_logs/%x_%A_%a.err" \
  --array="1-${BLAST_N}%6" scripts/05_blast_lsu_contigs.slurm
```

The reported contig analysis yielded results for 161/175
colonies overall and 146/158 CBASS colonies; these cohort counts are reported
analysis results and are not recalculated by the scripts in this checkout.

## Uniform read-level BLAST summary

[`scripts/06_compile_uniform_28s_summary.slurm`](scripts/06_compile_uniform_28s_summary.slurm)
is a single SLURM job that reads each sample's extracted R1 FASTQ, converts
FASTQ records to FASTA in a stream, and searches the LSU BLAST database with
BLAST+ `2.15.0+`, `-evalue 1e-10`, and `-max_target_seqs 5`. For each R1 read
with hits, it retains the highest-bitscore hit and counts *Durusdinium*
matches from that hit's title. It also adds the extracted paired-read count,
contig count, and best-scoring contig BLAST hit when available. Thus the
`pct_durusdinium_of_hits` column is the percentage among reads with a BLAST
hit, not among all reads. A row is emitted even for samples with no extracted
reads or no contigs.

Run it only after the clean-ID BLAST database, read-extraction outputs, and
contig BLAST tables are ready. The manifest is passed explicitly so the same
175-sample cohort is used:

```bash
mkdir -p "${SYMBIO_BASE_DIR}/slurm_logs"
SUMMARY_MANIFEST="${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt"
sbatch --chdir="${PWD}" \
  --export="ALL,SAMPLE_MANIFEST=${SUMMARY_MANIFEST}" \
  --output="${SYMBIO_BASE_DIR}/slurm_logs/%x_%j.out" \
  --error="${SYMBIO_BASE_DIR}/slurm_logs/%x_%j.err" \
  scripts/06_compile_uniform_28s_summary.slurm
```

The job writes
`blast_results/all_175_colonies_uniform_read_and_contig_calls.tsv` and copies
the summary and available LSU identifier map to
`~/symbiodiniaceae_lsu_backup/`. It does **not** delete `merged_fastp/`.
Inspect the SLURM log and table, validate the backup, and archive any other
required data before manually removing merged FASTQs. The reported analysis
included 12 low-coverage Tutuila CBASS colonies with 1-9 extracted read pairs
that fell below the assembly threshold; their reported best read hits were
98.5-100% identity to *Durusdinium*. These sequence matches are not
independent proof of species identity. The reported cohort totals (175
overall, 158 CBASS, including 146 with contigs) are analysis results, not
hard-coded assumptions in the job.

## Geneious FASTA export

[`scripts/07_export_geneious_fasta.sh`](scripts/07_export_geneious_fasta.sh)
combines every contig found for samples in the verified manifest with the 77
reference sequences. It renames colony contigs as
`>SAMPLE_ID|contig_<N>|len_<BP>` using the `len=` field in MEGAHIT headers.
Colonies without contigs are reported and skipped. The script refuses to
overwrite existing output or backup files; move specific prior exports before
rerunning.

Run the export after assemblies are complete, using the same cohort manifest:

```bash
EXPORT_MANIFEST="${SYMBIO_BASE_DIR}/sample_manifests/new_samples_only_manifest.txt"
srun -p standard -c 1 --mem=2G --time=00:10:00 \
  bash scripts/07_export_geneious_fasta.sh "${EXPORT_MANIFEST}"
```

It creates `geneious_export/all_161_colonies_28S_LSU_contigs.fasta` and
`geneious_export/all_161_colonies_plus_LaJeunesse2018_refs.fasta` under the
scratch base, and backs up both files to
`~/symbiodiniaceae_lsu_backup/geneious_export/`. The `all_161` filenames
reflect the reported cohort and do not enforce exactly 161 colonies or
contigs; the command prints actual sequence counts. The combined FASTA can be
used for MAFFT/Geneious alignment and tree inference. Geneious Prime `2026.1.3`
was reported for the downstream analysis; alignment parameters and tree
settings must also be archived to make the tree reproducible. The reference
panel spans seven genera but is not comprehensive or balanced, so species-
level interpretations need additional evidence.

## Storage, archiving, and publication acknowledgement

Sedna `/scratch` is unbacked and quota-limited. After validating all outputs,
archive the final contigs, extracted read counts, cohort manifests, BLAST
tables/database provenance, and Geneious FASTA/alignment/tree to
`/home/kturnham/symbiodiniaceae_lsu_backup/` (RAID6) and the project archive
at NOAA NCEI Accession `0318878`. Confirm the accession's current metadata and
contents before treating it as a complete backup. Delete `merged_fastp/`
FASTQs only after the pipeline no longer needs them and the final archived
inputs/results have been verified.

For publications, include the acknowledgement text currently required by
NOAA Fisheries for use of Sedna; use the official wording from the cluster
documentation rather than a paraphrase. The exact required wording is not
stored in this repository.
