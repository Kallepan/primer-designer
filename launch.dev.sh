#!/bin/bash

rm results/*
rm results/amplicons/*
rm tmp/*
rm logs/*
# python3 workflow/scripts/proto_primers.py -f data/myc_tuberculosis_h37Rv.fasta -p config/primer3_settings.yaml -r data/loci_formatted.csv -t tmp/ -o results/test.json
snakemake --use-conda -p results/myc_tuberculosis_h37Rv.proto_primers.json --cores 2