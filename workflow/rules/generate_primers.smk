import os

def get_amplicon_fasta_files():
    """ Returns a tuple of sample and amplicon names contained in the results/amplicons folder"""

    amplicons = []
    samples = []
    for root, dirs, files in os.walk("results/amplicons"):
        for file in files:
            if not file.endswith(".fasta"):
                continue
            
            samples.append(file.split(".")[0])
            amplicons.append(file.split(".")[1])
    
    return samples, amplicons

SAMPLES, AMPLICONS = get_amplicon_fasta_files() 

rule generate_proto_primers:
    """ Generates the proto-primers for each amplicon """
    input:
        expand(
            "results/amplicons/{sample}.{amplicon}.fasta",
            sample = SAMPLES,
            amplicon = AMPLICONS
        )
    conda:
        "envs/primers.yaml"
    output:
        "results/{sample}.proto"
    shell:
        "python3 scripts/generate_proto_primer.py {input}"
    