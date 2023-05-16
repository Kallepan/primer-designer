#!/bin/bash

# python3 workflow/scripts/filter_primers_by_alignment.py --alignment results/myc_tuberculosis_h37Rv.1.alignment.csv --primers results/myc_tuberculosis_h37Rv.1.proto_primers.json --output results/myc_tuberculosis_h37Rv.1.filtered_primers.json
# exit 0

rm results/*
rm results/* -rf
rm tmp/* -rf
rm logs/*

snakemake --use-conda -p results/myc_tuberculosis_h37Rv.summary.csv --cores 2