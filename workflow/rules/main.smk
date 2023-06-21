pool_count = config["primer_gen_config"]["pool_count"]
rule all:
    input: "logs/{species}.dump.log"
    output: temp("results/{species}.dummy")
    log: "results/{species}.dummy"
    shell:
        "echo 'dummy' > {output}"

rule dump_database:
    input: expand("results/{{species}}.{pool}.evaluated_primers.json", pool=range(0, pool_count))
    output: 
        "results/dump/{species}.proto_primers.tsv",
        "results/dump/{species}.alignments.tsv"
    log: "logs/{species}.dump.log"
    conda: "../envs/dump.yaml"
    shell:
        """
        python3 workflow/scripts/dump_database.py \
        --db results/{wildcards.species}.db \
        --species {wildcards.species} \
        --output_dir results/dump/ \
        &>> {log}"""