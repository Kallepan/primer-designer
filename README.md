# Primer Design Pipeline

## Requirements

### Config files

- Please use the config.yaml file and primer3_settings.yaml to adjust the settings of the pipeline

### Software requirements

- The pipeline requires the conda package manager, and snakemake to be installed
- Snakemake can be installed using conda. It will automatically install all other dependencies

## Description of Pipeline

### 1. Proto-Primers Generation

- Generate all possible primers using Primer3
- Allow the user to change primer3 config settings (tm, primer length, gc content)
- Simple wrapper written in python to call primer3
- Generation:
  - Two pools
  - Define regions of interest which are divided into amplicons
  - Each amplicons is flanked by forward and reverse primers
- Store output in json file
- Generate as many proto-primers as possible

#### Formula for determining if amplicons can overlap

Formula for calculating whether the amplicons can overlap or not.
If true then overlap is not possible:
$$A_{max} \leq 2 * A_{min} * (1 -min\_overlap)$$

#### PrimerGen Algorithm

- Generate seed sites with a certain distance between them
  - Distance between amplicon seed sites:
  $$min\_amp\_size*(1-min\_overlap*2) $$
  - pool offset: $$ min\_amp\_size*(1-min\_overlap) $$
- Generate primers for each amplicon using primer3:
  - Forward primer: 5'->3'
  - Reverse primer: 3'->5'
  - **primers are always stored and given in 5'->3' direction**
- Primers are generated for each amplicon
- Primers are only allowed in a so called buffer region:
  - Forward primers are placed in: 0 to buffer_region
  - Reverse primers are placed in: min_amplicon_size - buffer_region to min_amplicon_size
- If no forward or no reverse primers are found in this buffer region, the amplicon is skipped and not included in the final output
- Store the primers in the database with metadata extracted from the primer3 output

### 2. Filter Proto-Primers

#### Filter against target organism

- Align primers against the target organism using bowtie and store the results in a table
- Apply filters on the alignments:
  - The number of mismatches
  - The number of mismatches close to 3' end
  - The number of adjacent alignments of primers on the complementary strand
  - The number of times the primer misaligns to the genome
- Primers which do not meet filter criteria are marked as discarded in the database

#### Filter against other species

- Align primers against the target organism using bowtie and store the results in a table.
- Apply filters on the alignments:
  - The number of mismatches
  - The number of mismatches close to 3' end
  - The number of adjacent alignments of primers on the complementary strand
  - The number of times the primer misaligns to the genome
- All alignments are checked for adjacent alignments. The relationship of alignment with each other is stored in a graph. Using the minimal vertex cover problem, nodes/primers are selected which cover all edges/alignments. The selected primers are marked as discarded in the database.
- Bowtie uses pre-built indexes to align the primers against the genome. To reduce runtime, you can generate these indexes yourself and place them in the indexes folder. Alternatively genomes and indexes for bowtie are offered for download [Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn). Indexes can be generated using bowtie-build:
- Indexes should match in name with the given genome. If an index is not found, the pipeline will build the index using bowtie-build.

```bash
# How to build an index using bowtie-build:
bowtie-build foreign_species/sequence.fasta indexes/sequence
```

### 3. SADDLE

- Run the SADDLE algorithm on the filtered primers.
- SADDLE is a simulated annealing algorithm that optimizes primer sets for a given set of parameters.

#### SADDLE Algorithm

1. Generate a Set
2. Calculate loss for set (L = SUM(Badness(p1,p2))), where p1 and p2 are primers in the set
    - Badness is a function that calculates the badness of two primers
    - Badness takes into consideration the subsequences of the primers that could bind to other primers
    - Badness takes into consideration the distance of the subsequence to the 3' end of the primer
3. Generate a temporary set by replacing one or more primers from the current set by random replacement
4. Calculate loss for temporary set
5. Accept or reject temporary set based on a simulated annealing algorithm
6. Repeat 3-5 until the break off criteria is matched

```rust
if L(temporary_set) <= L(current_set) {
  current_set  = temporary_set
}
else {
  current_set = current_set     //(with probability 1-p)
  current_set = temporary_set   //(with probability p)
}
```

### 4. Report

- Display the results in a single html document (webbrowser) and export them in a .bed file (to be used with IGV) and a csv file

#### BED File

- List of primers with their positions in the genome

#### CSV

- List of amplicons with their respective forward and reverse primers in 5'->3' direction

#### HTML  

- Detailed view with sequences
- Simplified view with only an overview plot of the results
- Plot of loss of SADDLE function over time
- Overview of amplicons and the number of generated/discarded primers per amplicon in a table
- The html file is generated using angularjs, angular material, d3js, and gulp
