This repository contains pipelines and scripts supporting analyses with TelSPRITE, including the following:

-The directory called telsprite-pipeline contains files necessary for modifying the Guttman lab's SPRITE pipeline to create the TelSPRITE pipeline. Refer to the README file in that directory for specific instructions.

-The directory called tel-contact-frequencies has a small Snakemake pipeline that takes a cluster file with telomeres and extracts a plain text file containing telomere contacts in genomic bins. Refer to the README file in that directory for specific instructions.

-The directory called tel-sequences-ont contains a pipeline for calling ectopic telomere sequences from ONT long read WGS data. Refer to the README file in that directory for specific instructions.

-The script coverage_around_tel_breaks.py determines ONT sequencing coverage around ectopic telomere sequences, which can be outputted from the tel-sequences-ont pipeline.

-The script matrices_around_tel_breaks.py determines SPRITE sequencing coverage around ectopic telomere sequences.

-The script assembly_methods_mod.py is required to be present in the same directory to run coverage_around_tel_breaks.py and matrices_around_tel_breaks.py.

-The script normalized_and_cor_matrices.py takes SPRITE matrices for individual chromosomes as inputs and outputs Pearson correlation matrices of distance-normalized SPRITE contacts. These matrices are useful for visualizing A-B compartmentalization.

-The script tel_cluster_plot.py generates matrices with positions of telomere contacts from TelSPRITE data at high resolution, which can be used to generate cluster plots.
