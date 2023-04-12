import os

def get_amplicon_fasta_files():
    """ 
        Returns a tuple of species and amplicon names contained in the results/amplicons folder
        Currently not necessary, as the previous script generates a file with all amplicons
    """

    amplicons = []
    species = []
    for root, dirs, files in os.walk("results/amplicons"):
        for file in files:
            if not file.endswith(".fasta"):
                continue
            
            species.append(file.split(".")[0])
            amplicons.append(file.split(".")[1])
    
    return species, amplicons

SPECIES, AMPLICONS = get_amplicon_fasta_files() 

rule generate_proto_primers:
    """ Generates the proto-primers for each amplicon taking a fasta file with all amplicons as input"""
    input:
        "results/{species}.amplicons.fasta"
    conda:
        "../envs/primers.yaml"
    log:
        out = "logs/{species}.proto_primers.log"
    output:
        "results/{species}.proto_primers.json"
    shell:
        "python3 workflow/scripts/generate_proto_primer.py --input {input} --config config/primer3_settings.yaml > {log.out}"
