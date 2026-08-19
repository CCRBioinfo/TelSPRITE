# TelSPRITE: An approach for finding 3D telomere interactions with the genome

Most chromosome conformation studies omit repeat sequences that cannot be uniquely mapped. Because telomeres
become recombinogenic in cells that use Alternative Lengthening of Telomeres (ALT), we wished to understand
whether telomeres preferentially associate with certain regions of the genome, predisposing those regions to 
recombination with telomere sequence.  To achieve this end, we have made adaptations to the Guttman Lab SPRITE 
method to capture telomere sequences that are normally masked and include them in genome conformation 
analysis. 

For extensive information on installing and running the original Guttman Lab pipeline, see:

https://github.com/GuttmanLab/sprite-pipeline/wiki

Our tools layer onto this pipeline to provide additional telomere-related functionality.

This repository contains pipelines and scripts supporting analyses with TelSPRITE, including the following:

-The directory called telsprite-pipeline contains files necessary for modifying the Guttman lab's SPRITE pipeline to create the TelSPRITE pipeline. Refer to the README file in that directory for specific instructions.

-The directory called tel-contact-frequencies has a small Snakemake pipeline that takes a cluster file with telomeres and extracts a plain text file containing telomere contacts in genomic bins. Refer to the README file in that directory for specific instructions.

-The directory called tel-sequences-ont contains a pipeline for calling ectopic telomere sequences from ONT long read WGS data. Refer to the README file in that directory for specific instructions.

-The script coverage_around_tel_breaks.py determines ONT sequencing coverage around ectopic telomere sequences, which can be outputted from the tel-sequences-ont pipeline.

-The script matrices_around_tel_breaks.py determines SPRITE sequencing coverage around ectopic telomere sequences.

-The script assembly_methods_mod.py is required to be present in the same directory to run coverage_around_tel_breaks.py and matrices_around_tel_breaks.py.

-The script normalized_and_cor_matrices.py takes SPRITE matrices for individual chromosomes as inputs and outputs Pearson correlation matrices of distance-normalized SPRITE contacts. These matrices are useful for visualizing A-B compartmentalization.

-The script tel_cluster_plot.py generates matrices with positions of telomere contacts from TelSPRITE data at high resolution, which can be used to generate cluster plots.
