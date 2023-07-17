pools = config["metadata"]["pools"]
rule all:
    input:
        expand("{results}/dump/{{species}}.proto_primers.tsv", results=config["results_dir"]),
        expand("{results}/dump/{{species}}.alignments.tsv", results=config["results_dir"]),
        expand("{results}/{{species}}.bed", results=config["results_dir"]),
        expand("{results}/{{species}}.summary.html", results=config["results_dir"])
    output: touch("results/{species}.dummy")

chromosome = config["metadata"]["chromosome"]
rule format_into_bed:
    input: expand("{results}/{{species}}.primer_set.json", results=config["results_dir"])
    output: expand("{results}/{{species}}.bed", results=config["results_dir"])
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
    input: expand("results/{{species}}.{pool}.evaluated_primers.json", pool=pools)
    output: 
        expand("{results}/dump/{{species}}.proto_primers.tsv", results=config["results_dir"]),
        expand("{results}/dump/{{species}}.alignments.tsv", results=config["results_dir"])
    log: "logs/{species}.dump.log"
    conda: "../envs/base.yaml"
    params:
        output_dir = expand("{results}/dump", results=config["results_dir"])
    shell:
        """
        python3 workflow/scripts/dump_database.py \
            --db results/{wildcards.species}.db \
            --species {wildcards.species} \
            --output_dir {params.output_dir} \
        &> {log}"""