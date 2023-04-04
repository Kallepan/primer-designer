#!/bin/bash

rm output/*

python3 src/amplicons.py res/myc_tuberculosis_h37Rv.fasta res/loci_formatted.csv output/test.fasta
