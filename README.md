# Primer Design Pipeline

## 1. Proto-Primers Generation

- Generate all possible primers using Primer3
- Allow the user to change primer3 config settings (tm, primer length, gc content)
- Simple wrapper written in python to call primer3
- Generation:
  - Two pools
  - Define regions of interest which are divided into amplicons
  - Each amplicons is flanked by forward and reverse primers
- Store output in json file
- Generate as many proto-primers as possible

### Formula for determining if amplicons can overlap

```latex
\documentclass{article}
\begin{document}

Formula for calculating whether the amplicons can overlap or not.
If true then overlap is not possible:
\[ Amax <= 2*Amin*(1-min_overlap) \] 

\end{document}
```

### PrimerGen Algorithm

- Generate seed sites with a certain distance between them
  - Distance: min_amp_size * (1 - min_overlap)
  - pool offset: min_amp_size * (1 - min_overlap)
- Generate primers for each amplicon using primer3:
  - Forward primer: 5'->3'
  - Reverse primer: 3'->5' (output is in 5'->3' direction)
- Primers are generated for each amplicon
- Primers are only allowed in a so called buffer region:
  - Forward in: 0 to buffer_region
  - Reverse in: min_amplicon_size - buffer_region to min_amplicon_size
- If no primers are found in the buffer region -> skip amplicon
- Store the primers in the database with metadata

## 2. Filter Proto-Primers

### Filter against target organism

- Align primers against the target organism using bowtie
- Score primers
  - score each alignment based on the amount of mismatches
  - score each alignment based on the mismatches in the 3' end
  - score each alignment based on the amount of adjacent alignments
- Let the user choose the method of filtering. Either:
  - Hard filter -> discard all primers aligning multiple times with adjacent alignments
  - Soft filter -> score primers based on the amount of adjacent alignments
- Mark filtered primers in the database and store the score

### Filter against other species

- TODO: Align primers against other species using bowtie
- Simulate PCR? Taking into consideration mismatches and adjacent alignments
- To reduce runtime, you can generate these indexes yourself and place them in the indexes folder. Alternatively genomes and indexes are offered for download [Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn). Indexes can be generated using bowtie-build:

```bash
bowtie-build exclude/sequence.fasta indexes/sequence
```

## 3. SADDLE

- Run the SADDLE algorithm on the filtered primers.
- SADDLE is a simulated annealing algorithm that optimizes primer sets for a given set of parameters.

### SADDLE Algorithm

1. Generate a Set
2. Calculate loss for set (L = SUM(Badness(p1,p2))), where p1 and p2 are primers in the set
    - Badness is a function that calculates the badness of two primers
    - Badness takes into consideration the subsequences of the primers that could bind to other primers
    - Badness takes into consideration the distance of the subsequence to the 3' end of the primer
3. Generate a temporary set by replacing one or more primers from the current set by random replacement
4. Calculate loss for temporary set
5. Accept or reject temporary set based on a simulated annealing algorithm
6. Repeat 3-5 until convergence

```rust
if L(temporary_set) <= L(current_set) {
  current_set  = temporary_set
}
else {
  current_set = current_set     //(with probability 1-p)
  current_set = temporary_set   //(with probability p)
}
```

## 4. Report

- Display the results in a html document (webbrowser) and .bed file (IGV)
- Detailed view with sequences
- Simplified view with only an overview plot of the results
