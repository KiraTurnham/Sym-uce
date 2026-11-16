# Symbiodiniaceae LSU Read Extraction

This is a small Sedna SLURM workflow to enrich paired-end Illumina reads for
the Symbiodiniaceae large-subunit ribosomal RNA (LSU/28S) region, then assemble
the recovered pairs per sample. The initial goal is to check whether LSU
extraction and assembly work on holobiont data, including samples that may
contain multiple symbiont genera.

The pipeline does **not** make genus or species calls yet. It maps reads to a
multi-genus LSU reference, extracts properly paired reads, and assembles those
reads together per sample. Reference hits are enrichment baits only:
Bowtie2 reports one primary alignment for a read even when similar LSU
references in multiple genera are plausible. The downstream genus assessment
will be a separate step after reviewing the extracted-read and assembly
results. The current pipeline does not run BLAST.

## Requirements on Sedna

- SLURM and Environment Modules
- `aligners/bowtie2/2.5.0`
- `bio/samtools/1.19`
- Mamba/Conda `base` environment with MEGAHIT
- Bash, `gzip`, and standard GNU utilities

Confirm MEGAHIT is available before submitting jobs:

```bash
eval "$(mamba shell hook --shell bash)"
mamba activate base
megahit --version
```

## Setup

Clone the repository and create the scratch directories. Keep I/O-intensive
inputs, intermediate files, and outputs on `/scratch`, not NFS-backed home
storage.

```bash
git clone https://github.com/KiraTurnham/Sym-uce.git
cd Sym-uce

export SYMBIO_BASE_DIR="/scratch/${USER}/symbiodiniaceae_mining"
mkdir -p "${SYMBIO_BASE_DIR}"/{ref,sample_manifests,mapped_bam,assemblies,slurm_logs}
cp references/symbiodiniaceae_lsu.fasta \
  "${SYMBIO_BASE_DIR}/ref/symbiodiniaceae_lsu.fasta"
```

The reference contains 75 LSU sequences from the Symbiodiniaceae family in
the LaJeunesse et al. systematic revision alignment, plus two 28S records for
*Durusdinium*. Its seven represented genera include *Cladocopium*, so reads
from an unrepresented species such as a potential *Cladocopium* C40 can still
be recruited by genus-level LSU similarity. The reference does not establish
species identity. Provenance and limitations are documented in
[`references/README.md`](references/README.md) and
[`references/lsu_reference_sources.tsv`](references/lsu_reference_sources.tsv).

Set `FASTP_INPUT_DIR` to the existing `fastp_trimmed` directory. This workflow
expects paired files named `<sample>_R1_trimmed.fastq.gz` and
`<sample>_R2_trimmed.fastq.gz`, such as
`ESA_OCC_007_26_S11_L004_R1_trimmed.fastq.gz`. These suffixes are configured
in [`scripts/config.sh`](scripts/config.sh).

Build the Bowtie2 index on Sedna:

```bash
bash scripts/01_build_lsu_reference_index.sh
```

## Generic sample manifest

This manual option includes every matching sample in the directory. Do not use
it for a new-samples-only run if the directory also contains older samples.
Set `FASTP_INPUT_DIR` if it differs from the default, then produce a
one-sample-per-line manifest:

```bash
export FASTP_INPUT_DIR="/scratch/${USER}/phyluce_project/workflow_fast_megahit/fastp_trimmed"
find "${FASTP_INPUT_DIR}" -maxdepth 1 -type f -name '*_R1_trimmed.fastq.gz' \
  -printf '%f\n' | sed 's/_R1_trimmed\.fastq\.gz$//' \
  > "${SYMBIO_BASE_DIR}/sample_manifests/sample_list.txt"

cat "${SYMBIO_BASE_DIR}/sample_manifests/sample_list.txt"
```

Check that every listed sample has both mates before submission:

```bash
while IFS= read -r sample; do
  test -s "${FASTP_INPUT_DIR}/${sample}_R1_trimmed.fastq.gz" &&
    test -s "${FASTP_INPUT_DIR}/${sample}_R2_trimmed.fastq.gz" ||
    { echo "Missing FASTQ mate for ${sample}" >&2; exit 1; }
done < "${SYMBIO_BASE_DIR}/sample_manifests/sample_list.txt"
```

## Run

The default request is 8 CPUs, 16 GB RAM, and 3 hours per sample. The array is
throttled to six concurrent tasks. Ensure the log directory exists before
submitting.

```bash
NUM_SAMPLES=$(wc -l < "${SYMBIO_BASE_DIR}/sample_manifests/sample_list.txt")
if [[ "${NUM_SAMPLES}" -lt 1 ]]; then
  echo "Sample manifest is empty" >&2
  exit 1
fi

sbatch --output="${SYMBIO_BASE_DIR}/slurm_logs/%x_%A_%a.out" \
  --error="${SYMBIO_BASE_DIR}/slurm_logs/%x_%A_%a.err" \
  --array="1-${NUM_SAMPLES}%6" scripts/02_extract_and_assemble_lsu.slurm
```

The explicit log paths also work if `SYMBIO_BASE_DIR` is customized; SLURM
does not expand values sourced from `config.sh` inside `#SBATCH` directives.
Submit the job while your current directory is the repository root so the
script can source `scripts/config.sh` through SLURM's `SLURM_SUBMIT_DIR`.

## Pilot, then the remaining new samples

Make a manifest containing only the newly fastp-processed sample names,
without read suffixes:

```bash
cat > "${SYMBIO_BASE_DIR}/sample_manifests/new_samples.txt" <<'EOF'
NEW_SAMPLE_1
NEW_SAMPLE_2
EOF
```

Replace the example IDs with the actual sample names. Run one of them as a
pilot:

```bash
bash scripts/03_submit_lsu_one.sh NEW_SAMPLE_1
```

Review the pilot's log, recovered-read count, and assembly output before
submitting the rest. The all-samples submitter validates every sample name and
both non-empty FASTQ mates, snapshots the manifest under
`sample_manifests/`, and submits a SLURM array throttled to six concurrent
tasks. Its manifest must contain only samples that still need processing; if
the pilot succeeded, omit it from the batch list to avoid processing its
outputs again. The worker refuses to overwrite output files from any earlier
run.

```bash
grep -vx 'NEW_SAMPLE_1' \
  "${SYMBIO_BASE_DIR}/sample_manifests/new_samples.txt" \
  > "${SYMBIO_BASE_DIR}/sample_manifests/remaining_samples.txt"
bash scripts/04_submit_lsu_new_samples.sh \
  "${SYMBIO_BASE_DIR}/sample_manifests/remaining_samples.txt"
```

The submitter checks for the reference and a complete Bowtie2 index, creates
the log directory, and passes explicit log paths so customized scratch paths
work. The job requests 8 CPUs, 16 GB RAM, and 3 hours per sample. Set
`FASTP_INPUT_DIR` before running the submit helpers if it differs from the
default in `scripts/config.sh`.

## Outputs and interpretation

- `mapped_bam/<SAMPLE>/<SAMPLE>_lsu_R1.fastq.gz` and
  `<SAMPLE>_lsu_R2.fastq.gz`: properly paired reads aligned to one or more
  LSU references. Temporary BAM files are removed after extraction.
- `assemblies/<SAMPLE>/final.contigs.fa`: MEGAHIT contigs at least 300 bp.
- `slurm_logs/`: per-task SLURM standard output and error.

Samples with fewer than 10 recovered LSU read pairs are logged and skipped.
The recovered pairs are assembled together, not separated by genus. In a
mixed-genus sample, closely related templates may assemble separately or may
be collapsed; a contig is not by itself proof of a distinct genus or strain.
Review read counts, assembly contigs, and their alignments before adding a
separate BLAST-based genus-assessment step. Do not interpret a missing genus
hit as proof that genus is absent.
