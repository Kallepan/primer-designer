rule ensure_saddle_binary_is_build:
    input: "packages/saddle/cargo.toml"
    output: "packages/saddle/target/release/saddle"
    log: "logs/saddle/build.log"
    conda: "../envs/saddle.yaml"
    shell: "cd packages/saddle && cargo install && cargo build --release &> {log}"

rule run_saddle_for_each_pool:
    input: 
        binary = "packages/saddle/target/release/saddle",
        dump = "logs/{species}.dump.log",
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