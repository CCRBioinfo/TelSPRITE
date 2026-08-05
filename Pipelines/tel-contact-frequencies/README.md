This short pipeline extracts telomere contacts from TelSPRITE cluster files. The pipeline can be run as follows:

-Clone the pipeline

-Copy the scripts/HiCorrector_1.2 and scripts/python directories from the Guttman lab SPRITE pipeline (Quinodoz et al., 2022) into the scripts directory of this pipeline

-Replace scripts/python/assembly.py with the version of assembly.py found here

-Replace scripts/python/contact.py with the version of contact.py found here

-Ensure that you have generated cluster files from the telsprite-pipeline that contain telomere reads

-Edit config.yaml to include paths to appropriate files

-Run the pipeline by submitting the script run_pipeline.sh

Quinodoz, S.A., Bhat, P., Chovanec, P. et al. SPRITE: a genome-wide method for mapping higher-order 3D interactions in the nucleus using combinatorial split-and-pool barcoding. Nat Protoc 17, 36–75 (2022). https://doi.org/10.1038/s41596-021-00633-y