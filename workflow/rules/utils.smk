rule create_db:
    """ Creates a database file using sqlite3 """
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

rule export_amplicons:
    """ Exports information about each amplicon to a json file to be used for the report """
    input: 
        db = "results/{species}.db",
        regions = expand("results/{{species}}.{pool}.evaluated_primers.json", pool=pools),
    output:
        expand("{results}/{{species}}.amplicons.json", results = config["results_dir"])
    log: "logs/{species}.export_amplicons.log"
    conda: "../envs/base.yaml"
    shell:
        """
        python3 workflow/scripts/export_amplicons.py \
            --db {input.db} \
            --output {output} \
            &> {log}
        """

plotting_buffer = config["plot"]["plotting_buffer"]
rule export_regions:
    """ Exports the regions to a json file to be used for the report """
    input: 
        db = "results/{species}.db",
        regions = config["regions"],
        fasta = expand("{genomes_dir}/{target_genome}", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"])
    output: expand("{results}/{{species}}.regions.json", results = config["results_dir"])
    log: "logs/primer_gen/{species}.regions.log"
    conda: "../envs/biopython.yaml"
    params:
        plotting_buffer = plotting_buffer
    shell:
        """
        python3 workflow/scripts/export_regions.py \
            --db {input.db} \
            --regions {input.regions} \
            --fasta {input.fasta} \
            --plotting_buffer {params.plotting_buffer} \
            --output {output}  &> {log}
        """


rule format_primers_to_tsv:
    """ Formats the primers to a tsv file for easy inspection """
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



rule export_primers_to_bed:
    """ Exports the primers to a bed file for easy inspection in genome browsers such as IGV """
    input: 
        primer_set = expand("{results}/{{species}}.primer_set.json", results=config["results_dir"]),
        fasta = expand("{genomes_dir}/{target_genome}", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"])
    output: expand("{results}/{{species}}.bed", results=config["results_dir"])
    log: "logs/{species}.bed.log"
    conda: "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/export_primers_to_bed.py \
            --input {input.primer_set} \
            --fasta {input.fasta} \
            --output {output} \
        &> {log}"""


rule dump_database:
    """ Dumps all tables in the database to a tsv file for easy inspection """
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


rule export_evaluated_primers:
    """
        Exports the evaluated primers to a json file to be used by SADDLE (written in RUST).
        In the future SADDLE could adjusted to use the database directly
    """
    input:
        target_species = "results/filter/target/{species}.{pool}.alignment.eval.tsv",
        foreign_species = "results/filter/foreign/{species}.foreign_species.dummy",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.evaluated_primers.json"
    log: "logs/{species}.{pool}.export.log"
    conda: "../envs/primers.yaml"
    shell:
        """
        python3 workflow/scripts/export_evaluated_primers.py \
            --db {input.db} \
            --output {output} \
            --pool {wildcards.pool} &> {log}
        """