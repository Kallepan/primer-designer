import os

rule format_primers_into_fasta:
    input: "results/{species}.{pool}.proto_primers.json"
    output: "results/{species}.{pool}.proto_primers.fasta"
    log: "logs/{species}.{pool}.filtered_primers.log"
    conda:
        "../envs/primers.yaml"
    shell:
        "python workflow/scripts/format_into_fasta.py --input {input} --output {output} >> {log} 2>&1"

rule create_index_for_target:
    input: 
        fasta = "data/{species}.fasta"
    output:
        temp(expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 4))),
        temp(expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 2)))
    params:
        outdir = lambda w, input: os.path.join("tmp", "indexes", os.path.splitext(os.path.basename(input.fasta))[0])
    log: "logs/{species}.filtered_primers.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "bowtie-build {input.fasta} {params.outdir} >> {log} 2>&1"

rule align_primers_to_species:
    input:
        expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 4)),
        expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 2)),
        primers_fasta = "results/{species}.{pool}.proto_primers.fasta"
    output: "results/{species}.{pool}.alignment.csv"
    params:
        index = lambda w, input: os.path.join("tmp", "indexes", os.path.basename(input.primers_fasta).split(".")[0])
    log: "logs/{species}.{pool}.filtered_primers.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "python workflow/scripts/align_primers.py --primers {input.primers_fasta} --index {params.index} --output {output} >> {log} 2>&1"

rule filter_primers_by_alignment:
    input: 
        alignment = "results/{species}.{pool}.alignment.csv",
        original_primers = "results/{species}.{pool}.proto_primers.json"
    output: "results/{species}.{pool}.filtered_primers.json"
    log: "logs/{species}.{pool}.filtered_primers.log"
    conda:
        "../envs/primers.yaml"
    shell:
        "python workflow/scripts/filter_primers_by_alignment.py --alignment {input.alignment} --primers {input.original_primers} --output {output} >> {log} 2>&1"