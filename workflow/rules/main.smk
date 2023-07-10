pool_count = config["metadata"]["pool_count"]
rule all:
    input:
        "results/dump/{species}.proto_primers.tsv",
        "results/dump/{species}.alignments.tsv",
        "results/final/{species}.bed",
        "results/final/{species}.loss_set.json",
        "results/final/{species}.summary.html"
    output: temp("results/{species}.dummy")
    log: "logs/{species}.final.log"
    conda: "../envs/base.yaml"
    shell:
        "echo 'dummy' > {output}"

chromosome = config["metadata"]["chromosome"]
rule format_into_bed:
    input: "results/final/{species}.primer_set.json"
    output: "results/final/{species}.bed"
    log: "logs/{species}.bed.log"
    conda: "../envs/biopython.yaml"
    params: chromosome=chromosome
    shell:
        """
        python3 workflow/scripts/format_into_bed.py \
            --input {input} \
            --output {output} \
            --chrom {params.chromosome} \
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