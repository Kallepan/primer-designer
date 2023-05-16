import os

rule create_index_for_target:
    input:
        fasta = "data/{species}.fasta"
    output:
        temp(expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 5))),
        temp(expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 3)))
    params:
        outdir = lambda w, input: os.path.join("tmp", "indexes", os.path.splitext(os.path.basename(input.fasta))[0])
    log: "logs/{species}.index.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "bowtie-build {input.fasta} {params.outdir} >> {log} 2>&1"

rule all:
    input: expand("results/{{species}}.{pool}.proto_primers.fasta", pool=range(0, config["pool_count"]))
    output: "results/{species}.summary.csv"
    shell:
        "touch {output}"
        
rule format_pool_into_fasta:
    input: "results/{species}.db"
    output:
        file = "results/{species}.{pool}.proto_primers.fasta"
    log:
        file = "logs/{species}.{pool}.format.log"
    conda: "../envs/primers.yaml"
    shell:
        "python3 workflow/scripts/format_into_fasta.py --input {input} --output {output.file} --pool {wildcards.pool} &>> {log.file}"

rule align_primers_to_species:
    input:
        expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 5)),
        expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 3)),
        primers_fasta = "results/{species}.{pool}.proto_primers.fasta"
    output: "results/{species}.{pool}.alignment.csv"
    params:
        index = lambda w, input: os.path.join("tmp", "indexes", os.path.basename(input.primers_fasta).split(".")[0]),
        mismatches = config["mismatches"]
    log: "logs/{species}.{pool}.align.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "python3 workflow/scripts/align_primers.py --primers {input.primers_fasta} --index {params.index} --output {output} --mismatches {params.mismatches} &>> {log}"

rule filter_primers_by_alignment:
    input: 
        alignment = "results/{species}.{pool}.alignment.csv",
        original_primers = "results/{species}.{pool}.proto_primers.json"
    output: "results/{species}.{pool}.filtered_primers.json"
    log: "logs/{species}.{pool}.filter.log"
    params:
        adjacency_limit = config["adjacency_limit"]
        
    conda:
        "../envs/primers.yaml"
    shell:
        """ python3 workflow/scripts/filter_primers_by_alignment.py \
        --alignment {input.alignment} \
        --primers {input.original_primers} \
        --output {output} \
        --adjacency_limit {params.adjacency_limit} \
        &>> {log} """