#!/bin/bash
#
# This script is used to debug various parts of the pipeline.
# Its content should not be considered as part of the pipeline
# and should not be used as a reference for how to run the pipeline!

rm results/* -rf
rm logs/*
rm tmp/* -rf

snakemake --use-conda --use-singularity -p results/myc_tuberculosis_h37Rv.dummy -j 4
