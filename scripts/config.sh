#!/usr/bin/env bash

export BASE_SCRATCH="${SYMBIO_BASE_DIR:-/scratch/${USER:-$(id -un)}/symbiodiniaceae_mining}"
export FASTP_DIR="${FASTP_INPUT_DIR:-/scratch/${USER:-$(id -un)}/phyluce_project/workflow_fast_megahit/fastp_trimmed}"

export REF_DIR="${BASE_SCRATCH}/ref"
export REF_FASTA="${REF_DIR}/symbiodiniaceae_lsu.fasta"
export REF_INDEX_PREFIX="${REF_DIR}/sym_lsu_bt2"

export MANIFEST="${BASE_SCRATCH}/sample_manifests/sample_list.txt"
export BAM_DIR="${BASE_SCRATCH}/mapped_bam"
export ASSEMBLY_DIR="${BASE_SCRATCH}/assemblies"
export LOG_DIR="${BASE_SCRATCH}/slurm_logs"

export THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-8}}"
export MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-300}"
export MIN_MAPPED_PAIRS="${MIN_MAPPED_PAIRS:-10}"
