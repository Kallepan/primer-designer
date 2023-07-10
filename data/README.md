# Data folder

This folder contains the data used in the analysis. The following describes the files and their contents.

## Regions.csv

- name: name of the region (str)
- start: start of the region (int)
- end: end of the region (int)

E.G.:

```csv
name,start,end
chr1,0,1000000
chr2,1000000,2000000
```

## *.fasta

A fasta file containing the target sequence of the species you want the primers for.

## exclude/*.fasta

A collection of fasta files containing sequences that should be filtered against during primer design.
