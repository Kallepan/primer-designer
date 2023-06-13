rule create_tmp_dir:
    input: "data/{species}.fasta"
    log: "logs/{species}.proto_primers.log"
    conda: "../envs/primers.yaml"
    output: temp(directory("tmp/{species}/"))
    shell:
        "mkdir -p {output}"

pool_count = config["primer_gen_config"]["pool_count"]
amplicon_size = config["primer_gen_config"]["amplicon_size"]
amplicon_buffer = config["primer_gen_config"]["amplicon_buffer_size"]
rule generate_proto_primers:
    input:
        fasta = "data/{species}.fasta",
        primer3_config = "config/primer3_settings.yaml",
        tmp_dir = "tmp/{species}/",
        regions = "data/loci_formatted.csv",
        output_dir = "results/"
    params:
        amplicon_size = amplicon_size,
        amplicon_buffer = amplicon_buffer,
        pool_count = pool_count
    conda:
        "../envs/primers.yaml"
    output:
        "results/{species}.db"
    log:
        out = "logs/{species}.proto_primers.log"
    shell:
        """python3 workflow/scripts/proto_primers.py \
        -f {input.fasta} \
        -p {input.primer3_config} \
        -o {input.output_dir} \
        -r {input.regions} \
        -t {input.tmp_dir} >> {log.out} 2>&1"""