pool_count = config["metadata"]["pool_count"]
rule all:
    input:
        "results/dump/{species}.proto_primers.tsv",
        "results/dump/{species}.alignments.tsv",
        "results/{species}.bed",
        "results/{species}.primer_set.json"
    output: temp("results/{species}.dummy")
    log: "results/{species}.dummy"
    conda: "../envs/base.yaml"
    shell:
        "echo 'dummy' > {output}"

chromosome = config["metadata"]["chromosome"]
rule format_into_bed:
    input: "results/{species}.primer_set.json"
    output: "results/{species}.bed"
    log: "logs/{species}.bed.log"
    conda: "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/format_into_bed.py \
            --input {input} \
            --output {output} \
            --chrom chromosome \
        &> {log}"""

rule dump_database:
    input: expand("results/{{species}}.{pool}.evaluated_primers.json", pool=range(0, pool_count))
    output: 
        "results/dump/{species}.proto_primers.tsv",
        "results/dump/{species}.alignments.tsv"
    log: "logs/{species}.dump.log"
    conda: "../envs/base.yaml"
    shell:
        """
        python3 workflow/scripts/dump_database.py \
            --db results/{wildcards.species}.db \
            --species {wildcards.species} \
            --output_dir results/dump/ \
        &> {log}"""