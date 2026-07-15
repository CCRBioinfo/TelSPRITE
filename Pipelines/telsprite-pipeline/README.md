The original SPRITE pipeline is the work of the Guttman Lab. For information about the pipeline and instructions about how to clone it, refer to the following citation:

Quinodoz, S.A., Bhat, P., Chovanec, P. et al. SPRITE: a genome-wide method for mapping higher-order 3D interactions in the nucleus using combinatorial split-and-pool barcoding. Nat Protoc 17, 36–75 (2022). https://doi.org/10.1038/s41596-021-00633-y

After cloning the SPRITE pipeline, make the following modifications to create the TelSPRITE pipeline:

-Delete the Snakefile and replace with a copy of the Snakefile from this repository

-Delete config.yaml and replace with a copy of the config.yaml file from this repository. If necessary, change parameters such as telomere filtering cutoffs.

-Add a copy of run_telo.sh to the main directory of the pipeline. Edit the file to include a path to Bowtie2 indexes.

-Add a copy of filter_tel_reads.py to the subdirectory scripts/python.

To operate the pipeline with telomere filtering, run run_telo.sh.
