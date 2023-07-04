rule regions_to_json:
    input: 
        db = "results/{species}.db",
        regions = "data/regions.csv",
        fasta = "data/{species}.fasta"
    output: "results/{species}.regions.json"
    log: "logs/primer_gen/{species}.regions.log"
    conda: "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/regions_to_json.py \
            --db {input.db} \
            --regions {input.regions} \
            --fasta {input.fasta} \
            --output {output}  &> {log}
        """
        
pool_count = config["metadata"]["pool_count"]
amplicon_size = config["primer_gen_config"]["amplicon_size"]
amplicon_buffer = config["primer_gen_config"]["amplicon_buffer_size"]
rule generate_proto_primers:
    input:
        "results/{species}.regions.json",
        fasta = "data/{species}.fasta",
        primer3_config = "config/primer3_settings.yaml",
        db = "results/{species}.db"
    params:
        amplicon_size = amplicon_size,
        amplicon_buffer = amplicon_buffer,
        pool_count = pool_count,
        tmp_dir = temp("tmp/{species}/")
    conda: "../envs/primers.yaml"
    # Defining the log file as output is stupid, because it will be removed upon an error
    output: temp("results/{species}.proto_primers.dummy")
    log: "logs/primer_gen/{species}.proto_primers.log"
    shell:
        """
        python3 workflow/scripts/proto_primers.py \
            -f {input.fasta} \
            -p {input.primer3_config} \
            -d {input.db} \
            -t {params.tmp_dir} \
             &> {log} && touch {output}
        """