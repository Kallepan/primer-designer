plotting_buffer = config["plot"]["plotting_buffer"]
rule regions_to_json:
    input: 
        db = "results/{species}.db",
        regions = config["regions"],
        fasta = expand("{genomes_dir}/{target_genome}.fasta", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"])
    output: expand("{results}/{{species}}.regions.json", results = config["results_dir"])
    log: "logs/primer_gen/{species}.regions.log"
    conda: "../envs/biopython.yaml"
    params:
        plotting_buffer = plotting_buffer
    shell:
        """
        python3 workflow/scripts/regions_to_json.py \
            --db {input.db} \
            --regions {input.regions} \
            --fasta {input.fasta} \
            --plotting_buffer {params.plotting_buffer} \
            --output {output}  &> {log}
        """
        
pool_count = config["metadata"]["pool_count"]
min_overlap = config["primer_gen_config"]["min_overlap"]
min_amplicon_size = config["primer_gen_config"]["min_amplicon_size"]
max_amplicon_size = config["primer_gen_config"]["max_amplicon_size"]
rule generate_proto_primers:
    input:
        expand("{results}/{{species}}.regions.json", results = config["results_dir"]),
        fasta = expand("{genomes_dir}/{target_genome}.fasta", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"]),
        primer3_config = "config/primer3_settings.yaml",
        db = "results/{species}.db"
    params:
        pool_count = pool_count,
        min_overlap = min_overlap,
        min_amplicon_size = min_amplicon_size,
        max_amplicon_size = max_amplicon_size,
        tmp_dir = "tmp"
    conda: "../envs/primers.yaml"
    # Defining the log file as output is stupid, because it will be removed upon an error
    output: "results/{species}.proto_primers.dummy"
    log: "logs/primer_gen/{species}.proto_primers.log"
    shell:
        """
        python3 workflow/scripts/generate_proto_primers.py \
            -f {input.fasta} \
            -d {input.db} \
            -p {input.primer3_config} \
            -t {params.tmp_dir} \
            --min_overlap {params.min_overlap} \
            --min_amplicon_size {params.min_amplicon_size} \
            --max_amplicon_size {params.max_amplicon_size} \
            &> {log} && touch {output}
        """