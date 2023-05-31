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

rule all:
    input: expand("results/{{species}}.{pool}.evaluated_primers.json", pool=range(0, config["pool_count"]))
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
        mismatches = config["max_mismatches"]
    conda: "../envs/bowtie.yaml"
    shell:
        """
        python3 workflow/scripts/align_primers.py \
            --primers {input.primers_fasta} \
            --index {params.index} \
            --db {input.db} \
            --pool {wildcards.pool} \
            --mismatches {params.mismatches} \
            --output {output} &>> {log}
        """

rule score_alignments:
    input:
        alignment = "results/{species}.{pool}.alignment.csv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.alignment.scoring.csv"
    log: "logs/{species}.{pool}.alignment.scoring.log"
    conda: "../envs/primers.yaml"
    params:
        match_multiplier = config["match_multiplier"],
        primer_end_cutoff = config["primer_end_cutoff"],
        close_to_3_prime_factor = config["close_to_3_prime_factor"],
        far_from_3_prime_factor = config["far_from_3_prime_factor"],
        base_penalty = config["base_penalty"]
    shell:
        """python3 workflow/scripts/score_alignments.py \
        --db {input.db} \
        --output {output} \
        --pool {wildcards.pool} \
        --match_multiplier {params.match_multiplier} \
        --primer_end_cutoff {params.primer_end_cutoff} \
        --close_to_3_prime_factor {params.close_to_3_prime_factor} \
        --far_from_3_prime_factor {params.far_from_3_prime_factor} \
        --base_penalty {params.base_penalty} \
        &>> {log}"""

rule eval_primers_with_target:
    input:
        scoring = "results/{species}.{pool}.alignment.scoring.csv",
        db = "results/{species}.db",
    output: "results/{species}.{pool}.primer.scoring.csv"
    log: "logs/{species}.{pool}.evaluation.log"
    params:
        adjacency_limit = config["max_adjacency_limit"]
    conda:
        "../envs/primers.yaml"
    shell:
        """python3 workflow/scripts/eval_primers_with_target.py \
        --db {input.db} \
        --output {output} \
        --adjacency_limit {params.adjacency_limit} \
        --pool {wildcards.pool} &>> {log}"""

rule export_to_json:
    input:
        scoring = "results/{species}.{pool}.primer.scoring.csv",
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