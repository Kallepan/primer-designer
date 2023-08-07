import os


rule build_index:
    input: expand("{genomes_dir}/{{species}}.fasta", genomes_dir=config["genomes_dir"])
    output:
        expand("{index}/{{species}}.{version}.ebwt", version=range(1, 3), index=config["index_dir"]),
        expand("{index}/{{species}}.rev.{version}.ebwt", version=range(1, 3), index=config["index_dir"])
    params:
        outdir = lambda w, input: os.path.join(config["index_dir"], w.species)
    log: "logs/indexes/{species}.index.log"
    conda: "../envs/bowtie.yaml"
    shell: 
        """
        bowtie-build -r\
            {input} \
            {params.outdir} \
            &> {log}
        """


rule format_pool_into_fasta:
    input: 
        log = "logs/primer_gen/{species}.proto_primers.log",
        db = "results/{species}.db"
    output: "results/{species}.{pool}.proto_primers.fasta"
    log: "logs/filter/{species}.{pool}.format.log"
    conda: "../envs/primers.yaml"
    shell: 
        """
        python3 workflow/scripts/format_pool_into_fasta.py \
            --db {input.db} \
            --output {output} \
            --pool {wildcards.pool} \
            &> {log}
        """

max_mismatches = config["filter_settings"]["target_species"]["max_mismatches"]
rule align_primers_target_species:
    input:
        expand("{index}/{{species}}.{version}.ebwt", version=range(1, 3), index=config["index_dir"]),
        expand("{index}/{{species}}.rev.{version}.ebwt", version=range(1, 3), index=config["index_dir"]),
        primers_fasta = "results/{species}.{pool}.proto_primers.fasta",
        db = "results/{species}.db"
    output: "results/filter/target/{species}.{pool}.alignment.raw"
    log: "logs/filter/target/{species}.{pool}.alignment.log"
    params:
        index = lambda w, input: os.path.join(config["index_dir"], os.path.basename(input.primers_fasta).split(".")[0]),
        mismatches = max_mismatches
    conda: "../envs/bowtie.yaml"
    shell:
        """
        python3 workflow/scripts/align_primers.py \
            --primers {input.primers_fasta} \
            --index {params.index} \
            --mismatches {params.mismatches} \
            --output {output} &> {log}
        """


rule format_align_primers_target_species:
    input: 
        raw_alignment = "results/filter/target/{species}.{pool}.alignment.raw",
        db = "results/{species}.db"
    output: "results/filter/target/{species}.{pool}.alignment.tsv"
    log: "logs/filter/target/{species}.{pool}.alignment.format.log"
    conda: "../envs/biopython.yaml"
    shell: 
        """
        python3 workflow/scripts/format_align_primers.py \
            --input {input.raw_alignment} \
            --output {output} \
            --db {input.db} \
            --pool {wildcards.pool} \
            --species {wildcards.species} \
            &> {log}
        """


adjacency_limit = config["filter_settings"]["target_species"]["adjacency_limit"]
bases_to_ignore = config["filter_settings"]["target_species"]["bases_to_ignore"]
max_misalignments = config["filter_settings"]["target_species"]["max_misalignments"]
rule eval_primers_with_target_species:
    input:
        scores = "results/filter/target/{species}.{pool}.alignment.tsv",
        db = "results/{species}.db"
    output: "results/filter/target/{species}.{pool}.alignment.eval.tsv"
    log: "logs/filter/target/{species}.{pool}.alignment.eval.log"
    params:
        adjacency_limit = adjacency_limit,
        bases_to_ignore = bases_to_ignore,
        max_misalignments = max_misalignments
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 workflow/scripts/eval_primers_target.py \
            --db {input.db} \
            --output {output} \
            --pool {wildcards.pool} \
            --species {wildcards.species} \
            --adjacency_limit {params.adjacency_limit} \
            --bases_to_ignore {params.bases_to_ignore} \
            --max_misalignments {params.max_misalignments} \
            &> {log}
        """

