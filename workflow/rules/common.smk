""" This file contains all functions used by snakemake rules to run the pipeline """

import os


def get_indexes_input(w):
    """ Check if bowtie2 index is present, if not return bowtie index """
    if os.path.isfile(os.path.join(config["index_dir"], w.species + ".1.bt2")):
        return "bt2"
    else:
        return "ebwt"
        