# Primer Design Pipeline

## Description

- This pipeline is designed to design primers for a given genome and set of region of interests (ROIs) for Multiplex PCR reactions. It uses a tiled approach to design primers for a given set of ROIs. Each region is split into overlapping amplicons. Primers are designed for each amplicon and then filtered against the target genome and other species. The remaining primers are then optimized using the SADDLE algorithm.

## Requirements

### Config files

- Please use the config.yaml file and primer3_settings.yaml to adjust the settings of the pipeline
- For a detailed overview of how to run the pipeline please refer to: [.sample/README.md](.sample/README.md)

### Software

- [.sample/environment.yaml](.sample/environment.yaml) contains a list of all required software to run the pipeline. This pipeline uses the mamba environment manager to install all required software. If run using snakemake, the --use-conda flag will automatically install all required software using the environment.yaml file.
