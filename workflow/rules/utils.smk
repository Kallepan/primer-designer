rule create_db:
    input: 
        regions_sql_file = "workflow/sql/regions.sql",
        proto_primers_sql_file = "workflow/sql/proto_primers.sql",
        alignments_sql_file = "workflow/sql/alignments.sql",
    output: "results/{species}.db"
    log: "logs/{species}.db.log"
    conda: "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/create_db.py \
            --db {output} &> {log} \
            --regions_sql_file {input.regions_sql_file} \
            --proto_primers_sql_file {input.proto_primers_sql_file} \
            --alignments_sql_file {input.alignments_sql_file} \
            &> {log}
        """


plotting_buffer = config["plot"]["plotting_buffer"]
rule regions_to_json:
    input: 
        db = "results/{species}.db",
        regions = config["regions"],
        fasta = expand("{genomes_dir}/{target_genome}.fasta", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"])
    output: expand("{results}/{{species}}.regions.json", results = config["results_dir"])
    log: "logs/primer_gen/{species}.regions.log"
    conda: "../envs/biopython.yaml"
    params:
        plotting_buffer = plotting_buffer
    shell:
        """
        python3 workflow/scripts/regions_to_json.py \
            --db {input.db} \
            --regions {input.regions} \
            --fasta {input.fasta} \
            --plotting_buffer {params.plotting_buffer} \
            --output {output}  &> {log}
        """


rule format_primers_to_tsv:
    conda: "../envs/base.yaml"
    input: expand("{results}/{{species}}.primer_set.json", results=config["results_dir"])
    log: "logs/{species}.format_primers_to_tsv.log"
    output: expand("{results}/{{species}}.primer_set.tsv", results=config["results_dir"])
    shell:
        """
        python3 workflow/scripts/format_primers_to_tsv.py \
            --input {input} \
            --output {output} \
        &> {log}"""


chromosome = config["metadata"]["chromosome"]
rule export_primers_to_bed:
    input: expand("{results}/{{species}}.primer_set.json", results=config["results_dir"])
    output: expand("{results}/{{species}}.bed", results=config["results_dir"])
    log: "logs/{species}.bed.log"
    conda: "../envs/biopython.yaml"
    params: chromosome=chromosome
    shell:
        """
        python3 workflow/scripts/export_primers_to_bed.py \
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