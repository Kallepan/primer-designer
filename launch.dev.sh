#!/bin/bash

rm results/*
rm results/amplicons/*

snakemake --use-conda -p results/myc_tuberculosis_h37Rv.split --cores 1