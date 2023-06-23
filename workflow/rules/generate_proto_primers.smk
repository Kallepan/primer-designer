rule regions_to_json:
    input: 
        db = "results/{species}.db",
        regions = "data/regions.csv",
        fasta = "data/{species}.fasta"
    output: "results/{species}.regions.json"
    log: "logs/{species}.regions.log"
    conda: "../envs/converter.yaml"
    shell:
        """
        python3 workflow/scripts/regions_to_json.py \
            --db {input.db} \
            --regions {input.regions} \
            --fasta {input.fasta} \
            --output {output}  &>> {log}
        """
        
pool_count = config["primer_gen_config"]["pool_count"]
amplicon_size = config["primer_gen_config"]["amplicon_size"]
amplicon_buffer = config["primer_gen_config"]["amplicon_buffer_size"]
rule generate_proto_primers:
    input:
        fasta = "data/{species}.fasta",
        primer3_config = "config/primer3_settings.yaml",
        tmp_dir = "tmp/{species}/",
        regions = "results/{species}.regions.json",
        output_dir = "results/"
    params:
        amplicon_size = amplicon_size,
        amplicon_buffer = amplicon_buffer,
        pool_count = pool_count
    conda: "../envs/primers.yaml"
    output: "logs/{species}.proto_primers.log"
    log: "logs/{species}.proto_primers.log"
    shell:
        """python3 workflow/scripts/proto_primers.py \
        -f {input.fasta} \
        -p {input.primer3_config} \
        -o {input.output_dir} \
        -r {input.regions} \
        -t {input.tmp_dir} &>> {log}"""