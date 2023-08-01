# User Guide

- This guide briefly explains how to install and run the pipeline.

## Install conda environment

- Install conda environment from `environment.yml` file to install snakemake and other dependencies. This is necessary to ensure you have snakemake installed in the same environment as the pipeline. Snakemake will take care of installing all other dependencies in the environment locally. For more information on how to install conda, please refer to the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or the [mamba documentation](https://mamba.readthedocs.io/en/latest/installation.html).

```bash
# create the environment
mamba env create -f environment.yml

# activate the environment
mamba activate smk 
```

## Adjust the config files

- Adjust the config file `config/config.yaml` to adjust general settings of the pipeline. Please use the comments in the file to guide you.
- Adjust the `config/primer3_settings.yaml` file to adjust the settings of primer3 used to place primers. A detailed explanation of these configuration settings can be found [here](https://primer3.org/manual.html).

## Run the pipeline on a test dataset

- Run the pipeline on a test dataset to check if everything is working as expected.
- In the `.sample/data/` folder, there are two fasta files. `myc_tuberculosis_h37Rv.fasta` is the target organism for which we will design primers. `staph_aureus_nctc8325.fasta` is the organism against which our primers will be screened/filtered. Upon running to pipeline will automatically generate indexes for both fasta files.
- Note: for small genomes, the building of indexes is fast. For larger genomes such as human, this can take a while. If you are working with a large genome, you can download the bowtie indexes from [here](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn) and place them in the `.sample/indexes/` folder. The pipeline will automatically detect the indexes and use them instead of building them from scratch.

```bash
# run the pipeline
snakemake --use-conda --use-singularity -j 4 -p results/myc_tuberculosis_h37Rv.dummy
```

- Output files can be found in the `.sample/results/` folder.
- The `myc_tuberculosis_h37Rv.dummy` file is a dummy file that is used to trigger the pipeline. It is not used by the pipeline itself.
