""" This file contains all functions used by snakemake rules to run the pipeline """

import os


def get_align_primers_target_species_input(w):
    """ Check if bowtie2 index is present, if not return bowtie index """
    if os.path.isfile(os.path.join(config["index_dir"], w.species, ".1.bt2")):
        return [
            expand("{index}/{{species}}.{version}.bt2", version=range(1, 3), index=config["index_dir"]),
            expand("{index}/{{species}}.rev.{version}.bt2", version=range(1, 3), index=config["index_dir"]),
        ]
    else:
        return [
            expand("{index}/{{species}}.{version}.ebwt", version=range(1, 4), index=config["index_dir"]),
            expand("{index}/{{species}}.rev.{version}.ebwt", version=range(1, 4), index=config["index_dir"]),
        ]


def get_build_index_input(w):
    # Check if the .fasta file is present, if not return .fna
    if os.path.isfile(os.path.join(config["genomes_dir"], w.species, ".fasta")):
        return os.path.join(config["genomes_dir"], w.species, ".fasta")
    else:
        return os.path.join(config["genomes_dir"], w.species, ".fna")
