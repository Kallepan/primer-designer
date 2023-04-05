#!/bin/bash

rm results/*

python3 src/scripts/amplicons.py res/myc_tuberculosis_h37Rv.fasta res/loci_formatted.csv output/amplicons.fasta
