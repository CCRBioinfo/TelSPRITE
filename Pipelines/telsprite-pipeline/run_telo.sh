#!/bin/bash

source myconda
mamba activate base

module load snakemake || exit 1
module load samtools || exit 1
module load bowtie || exit 1
module load java || exit 1
module load cutadapt || exit 1
module load R || exit 1
module load picard || exit 1

export BOWTIE2_INDEXES=#Add path to Bowtie2 indexes

snakemake \
--snakefile Snakefile \
--use-conda \
-f \
-j 32 \
--cluster-config cluster.yaml \
--configfile config.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}" \
--rerun-incomplete \
telo

