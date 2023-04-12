#!/bin/bash

rm results/*
rm results/amplicons/*
rm tmp/*

snakemake --use-conda -p results/myc_tuberculosis_h37Rv.proto_primers.json --cores 1