#!/bin/bash

rm results/*

snakemake --use-conda -p results/myc_tuberculosis_h37Rv.amplicons.fasta --cores 1