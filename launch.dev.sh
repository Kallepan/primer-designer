#!/bin/bash
#
# This script is used to debug various parts of the pipeline.
# Its content should not be considered as part of the pipeline
# and should not be used as a reference for how to run the pipeline!

rm /homes/kkandeepan/data/results/* -rf
rm results -rf
rm logs -rf
rm tmp -rf

# Snakemake: use conda, use singularity, 4 cores, print shell commands
snakemake --use-conda --use-singularity -j 4 -p results/myc_tuberculosis_h37Rv.dummy