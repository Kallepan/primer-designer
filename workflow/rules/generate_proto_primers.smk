min_overlap = config["primer_gen_config"]["min_overlap"]
min_amplicon_size = config["primer_gen_config"]["min_amplicon_size"]
max_amplicon_size = config["primer_gen_config"]["max_amplicon_size"]
pool_count = config["metadata"]["pool_count"]
rule generate_proto_primers:
    input:
        expand("{results}/{{species}}.regions.json", results = config["results_dir"]),
        fasta = expand("{genomes_dir}/{target_genome}", genomes_dir = config["genomes_dir"], target_genome = config["target_genome"]),
        primer3_config = "config/primer3_settings.yaml",
        db = "results/{species}.db",
    params:
        min_overlap = min_overlap,
        min_amplicon_size = min_amplicon_size,
        max_amplicon_size = max_amplicon_size,
        pool_count = pool_count,
    conda: "../envs/primers.yaml"
    output: 
        touch("results/{species}.proto_primers.dummy"),
        output_dir = directory("results/{species}/amplicons/"),
    log: "logs/primer_gen/{species}.proto_primers.log"
    shell:
        """
        python3 workflow/scripts/generate_proto_primers.py \
            -f {input.fasta} \
            -d {input.db} \
            -p {input.primer3_config} \
            -o {output.output_dir} \
            --min_overlap {params.min_overlap} \
            --min_amplicon_size {params.min_amplicon_size} \
            --max_amplicon_size {params.max_amplicon_size} \
            --pool_count {params.pool_count} \
            &> {log}
        """