import os

rule create_index_for_target:
    input:
        fasta = "data/{species}.fasta"
    output:
        temp(expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 5))),
        temp(expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 3)))
    params:
        outdir = lambda w, input: os.path.join("tmp", "indexes", os.path.splitext(os.path.basename(input.fasta))[0])
    log: "logs/{species}.index.log"
    conda: "../envs/bowtie.yaml"
    shell:
        "bowtie-build {input.fasta} {params.outdir} >> {log} 2>&1"

pool_count = config["primer_gen_config"]["pool_count"]
rule all:
    input: expand("results/{{species}}.{pool}.evaluated_primers.json", pool=range(0, pool_count))
    output: "results/{species}.summary.csv"
    shell:
        "touch {output}"
        
rule format_pool_into_fasta:
    input: "results/{species}.db"
    output:
        file = "results/{species}.{pool}.proto_primers.fasta"
    log:
        file = "logs/{species}.{pool}.format.log"
    conda: "../envs/primers.yaml"
    shell:
        "python3 workflow/scripts/format_into_fasta.py --db {input} --output {output.file} --pool {wildcards.pool} &>> {log.file}"

max_mismatches = config["alignment_settings"]["max_mismatches"]
rule align_primers_to_species:
    input:
        expand("tmp/indexes/{{species}}.{version}.ebwt", version=range(1, 5)),
        expand("tmp/indexes/{{species}}.rev.{version}.ebwt", version=range(1, 3)),
        primers_fasta = "results/{species}.{pool}.proto_primers.fasta",
        db = "results/{species}.db"
    output: "results/{species}.{pool}.alignment.csv"
    log: "logs/{species}.{pool}.alignment.log"
    params:
        index = lambda w, input: os.path.join("tmp", "indexes", os.path.basename(input.primers_fasta).split(".")[0]),
        mismatches = max_mismatches
    conda: "../envs/bowtie.yaml"
    shell:
        """
        python3 workflow/scripts/align_primers.py \
            --primers {input.primers_fasta} \
            --index {params.index} \
            --db {input.db} \
            --pool {wildcards.pool} \
            --species {wildcards.species} \
            --mismatches {params.mismatches} \
            --output {output} &>> {log}
        """

base_penalty = config["alignment_settings"]["base_penalty"]
alignment_weight = config["alignment_settings"]["alignment_weight"]
mismatch_weight = config["alignment_settings"]["mismatch_weight"]
rule score_alignments:
    input:
        alignment = "results/{species}.{pool}.alignment.csv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.alignment.scores.csv"
    log: "logs/{species}.{pool}.alignment.scores.log"
    conda: "../envs/primers.yaml"
    params:
        base_penalty = base_penalty,
        alignment_weight = alignment_weight,
        mismatch_weight = mismatch_weight,
    shell:
        """python3 workflow/scripts/score_alignments.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} \
        --species {wildcards.species} \
        --alignment_weight {params.alignment_weight} \
        --mismatch_weight {params.mismatch_weight} \
        --base_penalty {params.base_penalty} \
        &>> {log}"""

adjacency_limit = config["evaluation_settings"]["max_adjacency_limit"]
hard_filter = config["evaluation_settings"]["hard_filter"]
rule eval_primers_with_target:
    input:
        scores = "results/{species}.{pool}.alignment.scores.csv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.proto_primers.scores.csv"
    log: "logs/{species}.{pool}.evaluation.log"
    params:
        adjacency_limit = adjacency_limit,
        hard_filter = hard_filter,
    conda:
        "../envs/primers.yaml"
    shell:
        # Check if we are using hard or soft filtering
        """
        if [ "{params.hard_filter}" == "True" ]; then
            python3 workflow/scripts/eval_primers_hard.py \
            --db {input.db} \
            --output {output} \
            --pool {wildcards.pool} \
            --species {wildcards.species} \
            --adjacency_limit {params.adjacency_limit} \
            &>> {log}
        else
            python3 workflow/scripts/eval_primers_soft.py \
            --db {input.db} \
            --output {output} \
            --pool {wildcards.pool} \
            --species {wildcards.species} \
            --adjacency_limit {params.adjacency_limit} \
            &>> {log}
        fi
        """

rule export_to_json:
    input:
        scores = "results/{species}.{pool}.proto_primers.scores.csv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.evaluated_primers.json"
    log: "logs/{species}.{pool}.export.log"
    conda:
        "../envs/primers.yaml"
    shell:
        """python3 workflow/scripts/format_into_json.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} &>> {log}"""