# User Guide

- This guide briefly explains how to install and run the pipeline.

## Install conda environment

- Install the mamba package manager. Mamba is a drop-in replacement for conda. It is not necessary to use mamba, but it is recommended as it is much faster than conda.
- Create a conda environment called `smk` and activate it using the commands below.
- This is necessary to ensure you have snakemake installed. Snakemake will take care of installing all other dependencies in the environment locally. For more information on how to install conda, please refer to the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or the [mamba documentation](https://mamba.readthedocs.io/en/latest/installation.html).

```bash
# create the environment (preferred)
mamba create -c conda-forge -c bioconda --name smk snakemake
# or run
mamba env create -f .sample/environment.yml

# activate the environment
mamba activate smk 
```

## Adjust the config files

- Adjust the config file `config/config.yaml` to adjust general settings of the pipeline. Please use the comments in the file to guide you.
- Adjust the `config/primer3_settings.yaml` file to adjust the settings of primer3 used to place primers. A detailed explanation of these configuration settings can be found [here](https://primer3.org/manual.html).

## Run the pipeline on a test dataset

- Run the pipeline on a test dataset to check if everything is working as expected.
- In the `.sample/data/` folder, there are two fasta files. `myc_tuberculosis_h37Rv.fasta` [(source)](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3) is the target organism for which we will design primers. `staph_aureus_nctc8325.fasta` [(source)](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1280/) is the organism against which our primers will be screened/filtered. Fasta files of genomes can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov). Upon running to pipeline will automatically generate indexes for both fasta files.
- Note: for small genomes, the building of indexes is fast. For larger genomes such as human, this can take a while. If you are working with a large genome, you can download the bowtie indexes from [here](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn) and place them in the `.sample/indexes/` folder. If you, place the indexes into the indexes folder, the beginning of the filename should match the beginning of the filename of the fasta files:

```bash
# example
data/GRCh38.fasta
indexes/GRCh38.1.ebwt
indexes/GRCh38.2.ebwt
indexes/GRCh38.rev.1.ebwt
indexes/GRCh38.rev.2.ebwt
```

The pipeline will automatically detect the indexes and use them instead of building them from scratch.

```bash
# run the pipeline
snakemake --use-conda --use-singularity -j 4 -p results/myc_tuberculosis_h37Rv.dummy
```

- Output files can be found in the `.sample/results/` folder.
- The `myc_tuberculosis_h37Rv.dummy` file is a dummy file that is used to trigger the pipeline. It is not used by the pipeline itself.
