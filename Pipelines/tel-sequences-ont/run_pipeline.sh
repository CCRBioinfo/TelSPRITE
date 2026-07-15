#!/bin/bash

module load snakemake
module load python
#module load minimap2
#module load samtools

snakemake --cores --configfile config.yaml

