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
        regions = "data/loci_formatted.csv",
        output_dir = "results/"
    params:
        amplicon_size = config["amplicon_size"],
        amplicon_buffer = config["amplicon_buffer_size"],
        pool_count = config["pool_count"]
    conda:
        "../envs/primers.yaml"
    output: 
        files = expand("results/{{species}}.{pool}.proto_primers.json", pool=range(1, config["pool_count"] + 1)),
        file = "results/{species}.proto_primers.json"
    log:
        out = "logs/{species}.proto_primers.log"
    shell:
        "python3 workflow/scripts/proto_primers.py -f {input.fasta} -p {input.primer3_config} -o {input.output_dir} -r {input.regions} -t {input.tmp_dir} >> {log.out} 2>&1"