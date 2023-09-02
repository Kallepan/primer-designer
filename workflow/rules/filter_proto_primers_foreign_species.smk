fasta_files = [os.path.join(config["genomes_dir"], genome) for genome in config["foreign_species_fastas"]]
fasta_names = [os.path.splitext(os.path.basename(f))[0] for f in fasta_files]


rule all_foreign_species:
    input:
        expand("results/filter/foreign/{{species}}.{fasta}.{pool}.alignment.eval.tsv", fasta=fasta_names, pool=pools)
    output: touch("results/filter/foreign/{species}.foreign_species.dummy")


max_mismatches = config["filter_settings"]["foreign_species"]["max_mismatches"]
rule align_primers_foreign_species:
    input:
        # ancient is a helper function that ensures the index is not rebuilt if it already exists
        lambda w: ancient(expand("{index}/{{species}}.{version}.{ending}", version=range(1, 3), index=config["index_dir"], ending=get_indexes_input(w))),
        lambda w: ancient(expand("{index}/{{species}}.rev.{version}.{ending}", version=range(1, 3), index=config["index_dir"], ending=get_indexes_input(w))),
        primers_fasta="results/{species}.{pool}.proto_primers.fasta",
        db = "results/{species}.db"
    output: "results/filter/foreign/{species}.{fasta}.{pool}.alignment.raw"
    log: "logs/filter/foreign/{species}.{fasta}.{pool}.alignment.log"
    params:
        index = lambda w: os.path.join(config["index_dir"], w.fasta),
        mismatches = max_mismatches,
    conda: "../envs/bowtie.yaml"
    shell:
        """
        python3 workflow/scripts/align_primers.py \
            --primers {input.primers_fasta} \
            --index {params.index} \
            --mismatches {params.mismatches} \
            --output {output} &> {log}
        """


rule format_align_primers_foreign_species:
    input: 
        raw_alignment = "results/filter/foreign/{species}.{fasta}.{pool}.alignment.raw",
        db = "results/{species}.db"
    output: "results/filter/foreign/{species}.{fasta}.{pool}.alignment.tsv"
    log: "logs/filter/foreign/{species}.{fasta}.{pool}.alignment.format.log"
    conda: "../envs/biopython.yaml"
    shell: 
        """
        python3 workflow/scripts/format_align_primers.py \
            --input {input.raw_alignment} \
            --output {output} \
            --db {input.db} \
            --pool {wildcards.pool} \
            --species {wildcards.fasta} \
            &> {log}
        """


adjacency_limit = config["filter_settings"]["foreign_species"]["adjacency_limit"]
bases_to_ignore = config["filter_settings"]["foreign_species"]["bases_to_ignore"]
rule eval_align_primers_foreign_species:
    input: 
        alignment = "results/filter/foreign/{species}.{fasta}.{pool}.alignment.tsv",
        db = "results/{species}.db"
    output: "results/filter/foreign/{species}.{fasta}.{pool}.alignment.eval.tsv"
    log: "logs/filter/foreign/{species}.{fasta}.{pool}.alignment.eval.log"
    conda: "../envs/base.yaml"
    threads: workflow.cores * 0.5
    params:
        adjacency_limit = adjacency_limit,
        bases_to_ignore = bases_to_ignore
    shell:
        """
            python3 workflow/scripts/eval_primers_foreign.py \
            --db {input.db} \
            --pool {wildcards.pool} \
            --species {wildcards.fasta} \
            --output {output} \
            --adjacency_limit {params.adjacency_limit} \
            --bases_to_ignore {params.bases_to_ignore} \
            &> {log}
        """

