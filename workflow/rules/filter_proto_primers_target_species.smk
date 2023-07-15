import os

rule build_index:
    input: expand("{genomes_dir}/{{species}}.fasta", genomes_dir=config["genomes_dir"])
    output:
        expand("{index}/{{species}}.{version}.ebwt", version=range(1, 5), index=config["index_dir"]),
        expand("{index}/{{species}}.rev.{version}.ebwt", version=range(1, 3), index=config["index_dir"])
    params:
        outdir = lambda w, input: os.path.join(config["index_dir"], w.species)
    log: "logs/indexes/{species}.index.log"
    conda: "../envs/bowtie.yaml"
    shell: "bowtie-build {input} {params.outdir} &> {log}"

rule format_pool_into_fasta:
    input: 
        log = "logs/primer_gen/{species}.proto_primers.log",
        db = "results/{species}.db"
    output: "results/{species}.{pool}.proto_primers.fasta"
    log: "logs/filter/{species}.{pool}.format.log"
    conda: "../envs/primers.yaml"
    shell: "python3 workflow/scripts/format_into_fasta.py --db {input.db} --output {output} --pool {wildcards.pool} &> {log}"

max_mismatches = config["alignment_settings"]["max_mismatches"]
rule align_primers_to_species:
    input:
        expand("{index}/{{species}}.{version}.ebwt", version=range(1, 5), index=config["index_dir"]),
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

rule format_align_primers_to_species:
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

base_penalty = config["alignment_settings"]["base_penalty"]
alignment_weight = config["alignment_settings"]["alignment_weight"]
mismatch_weight = config["alignment_settings"]["mismatch_weight"]
rule score_alignments:
    input:
        alignment = "results/filter/target/{species}.{pool}.alignment.tsv",
        db = "results/{species}.db",
    output: "results/filter/target/{species}.{pool}.alignment.scores.tsv"
    log: "logs/filter/target/{species}.{pool}.alignment.scores.log"
    conda: "../envs/primers.yaml"
    params:
        base_penalty = base_penalty,
        alignment_weight = alignment_weight,
        mismatch_weight = mismatch_weight,
    shell: """
        python3 workflow/scripts/score_alignments.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} \
        --species {wildcards.species} \
        --alignment_weight {params.alignment_weight} \
        --mismatch_weight {params.mismatch_weight} \
        --base_penalty {params.base_penalty} \
        &> {log}
    """

adjacency_limit = config["evaluation_settings"]["target_filter"]["max_adjacency_limit"]
badness_threshold = config["evaluation_settings"]["target_filter"]["threshold"]
rule eval_primers_with_target:
    input:
        scores = "results/filter/target/{species}.{pool}.alignment.scores.tsv",
        db = "results/{species}.db",
    output: "results/filter/target/{species}.{pool}.proto_primers.scores.tsv"
    log: "logs/filter/target/{species}.{pool}.evaluation.log"
    params:
        adjacency_limit = adjacency_limit,
        badness_threshold = badness_threshold
    conda:
        "../envs/primers.yaml"
    shell:
        """
        python3 workflow/scripts/eval_primers_target.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} \
        --species {wildcards.species} \
        --adjacency_limit {params.adjacency_limit} \
        --badness_threshold {params.badness_threshold} \
        &> {log}
        """

rule export_to_json:
    input:
        scores = "results/filter/target/{species}.{pool}.proto_primers.scores.tsv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.evaluated_primers.json"
    log: "logs/{species}.{pool}.export.log"
    conda: "../envs/primers.yaml"
    shell:
        """
        python3 workflow/scripts/format_into_json.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} &> {log}
        """
