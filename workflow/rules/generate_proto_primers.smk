rule create_tmp_dir:
    input:
        "data/{species}.fasta"
    log:
        out = "logs/{species}.proto_primers.log"
    output:
        temp(directory("tmp/{species}/"))
    shell:
        "mkdir -p {output}"

rule generate_proto_primers:
    input:
        fasta = "data/{species}.fasta",
        primer3_config = "config/primer3_settings.yaml",
        tmp_dir = "tmp/{species}/",
        regions = "data/loci_formatted.csv"
    params:
        amplicon_size = config["amplicon_size"],
        amplicon_buffer = config["amplicon_buffer_size"]
    conda:
        "../envs/primers.yaml"
    output:
        file = "results/{species}.proto_primers.json"
    log:
        out = "logs/{species}.proto_primers.log"
    shell:
        "python3 workflow/scripts/proto_primers.py -f {input.fasta} -p {input.primer3_config} -o {output.file} -r {input.regions} -t {input.tmp_dir} > {log.out} 2>&1"