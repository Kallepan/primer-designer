rule ensure_saddle_binary_is_build:
    input: "packages/saddle/cargo.toml"
    output: "packages/saddle/target/release/saddle"
    log: "logs/saddle/build.log"
    conda: "../envs/saddle.yaml"
    shell: "cd packages/saddle && cargo install && cargo build --release &> {log}"

rule run_saddle:
    input: 
        binary = "packages/saddle/target/release/saddle",
        pool = "results/{species}.{pool}.evaluated_primers.json"
    output: 
        set = "results/{species}.{pool}.saddle.set.json",
        loss = "results/{species}.{pool}.saddle.loss.json"
    log: "logs/saddle/{species}.{pool}.log"
    conda: "../envs/saddle.yaml"
    shell:
        """
        {input.binary} \
            --input-file {input.pool} \
            --output-file-set {output.set} \
            --output-file-loss {output.loss} \
            &> {log}
        """

pool_count = config["primer_gen_config"]["pool_count"]
species = config["primer_gen_config"]["species"]
rule merge_saddle_output:
    input: expand("results/{species}.{pool}.saddle.set.json", species=species,pool=range(0, pool_count))
    conda: "../envs/base.yaml"
    log: "logs/saddle/merge.log"
    output: "results/saddle.set.json"
    shell:
        """
        python workflow/scripts/merge_saddle_output.py \
            --output {output} \
            --input {input} \
            &> {log}
        """