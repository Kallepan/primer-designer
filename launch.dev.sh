#!/bin/bash
#
# This script is used to debug various parts of the pipeline.
# Its content should not be considered as part of the pipeline
# and should not be used as a reference for how to run the pipeline!

rm results/* -rf
rm logs/*

snakemake --use-conda -p results/myc_tuberculosis_h37Rv.summary.csv --cores 2

cd helper_scripts
python3 dump.py
