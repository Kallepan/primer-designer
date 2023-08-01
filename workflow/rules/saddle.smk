rule build_saddle_binary:
    input: "packages/saddle/Cargo.toml"
    output: "packages/saddle/target/release/saddle"
    log: "logs/saddle/build.log"
    conda: "../envs/saddle.yaml"
    shell: "cargo build --manifest-path {input} --release &> {log}"


min_subsequence_size = config["saddle"]["min_subsequence_size"]
max_subsequence_size = config["saddle"]["max_subsequence_size"]
max_iterations = config["saddle"]["max_iterations"]
amplicons_weight = config["saddle"]["amplicons_weight"]
primers_weight = config["saddle"]["primers_weight"]
rule run_saddle:
    input: 
        binary = "packages/saddle/target/release/saddle",
        pool = "results/{species}.{pool}.evaluated_primers.json"
    output:
        set = "results/{species}.{pool}.saddle_set.json",
        loss = "results/{species}.{pool}.saddle_loss.json"
    params:
        min_subsequence_size = min_subsequence_size,
        max_subsequence_size = max_subsequence_size,
        max_iterations = max_iterations,
        amplicons_weight = amplicons_weight,
        primers_weight = primers_weight
    log: "logs/saddle/{species}.{pool}.log"
    conda: "../envs/saddle.yaml"
    shell:
        """
        {input.binary} \
            --input-file {input.pool} \
            --output-file-set {output.set} \
            --output-file-loss {output.loss} \
            --min-subsequence-size {params.min_subsequence_size} \
            --max-subsequence-size {params.max_subsequence_size} \
            --max-iterations {params.max_iterations} \
            --amplicons-weight {params.amplicons_weight} \
            --primers-weight {params.primers_weight} \
            &> {log}
        """


rule merge_saddle_data_output:
    input: expand("results/{{species}}.{pool}.saddle_set.json", pool=pools)
    conda: "../envs/base.yaml"
    log: "logs/saddle/{species}.merge.log"
    output: expand("{results}/{{species}}.primer_set.json", results=config["results_dir"])
    shell:
        """
        python workflow/scripts/merge_saddle_output.py \
            --output {output} \
            --input {input} \
            &> {log}
        """


rule merge_saddle_loss_output:
    input: expand("results/{{species}}.{pool}.saddle_loss.json", pool=pools)
    conda: "../envs/base.yaml"
    log: "logs/saddle/{species}.merge.log"
    output: expand("{results}/{{species}}.loss_set.json", results=config["results_dir"])
    shell:
        """
        python workflow/scripts/merge_saddle_output.py \
            --output {output} \
            --input {input} \
            &> {log}
        """

