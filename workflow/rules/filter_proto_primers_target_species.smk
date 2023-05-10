
rule format_primers_into_fasta:
    input: "results/{species}.{pool}.proto_primers.json"
    output: "results/{species}.{pool}.proto_primers.fasta"
    log: "logs/{species}.{pool}.filtered_primers.log"
    conda:
        "../envs/primers.yaml"
    shell:
        "python workflow/scripts/format_into_fasta.py --input {input} --output {output} >> {log} 2>&1"

rule create_index_for_target_species:
    input: "data/{species}.fasta"
    output:
        temp(expand("tmp/index/{{species}}.{version}.ebwt", version=range(1, 4))),
        temp(expand("tmp/index/{{species}}.rev.{version}.ebwt", version=range(1, 2)))
    params:
        outdir = "tmp/index/"
    log: "logs/{species}.filtered_primers.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "bowtie-build {input} {params.outdir}/{wildcards.species} >> {log} 2>&1"

rule align_primers_to_species:
    input:
        expand("tmp/index/{{species}}.{version}.ebwt", version=range(1, 4)),
        expand("tmp/index/{{species}}.rev.{version}.ebwt", version=range(1, 2)),
        primers_fasta = "results/{species}.{pool}.proto_primers.fasta"
    output: "results/{species}.{pool}.filtered_primers.json"
    params:
        index = "tmp/index/{species}"
    log: "logs/{species}.{pool}.filtered_primers.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "python scripts/align_primers.py --primers_input {input.primers_fasta} --index_input {params.index} --output_file {output} >> {log} 2>&1"