# Fetch all fasta files to be excluded against
fasta_files = config["foreign_species_genomes"]
foreign_species = [os.path.splitext(os.path.basename(f))[0] for f in fasta_files]

rule all_foreign_species:
    # TODO: This rule is broken. It should be a union of all the eval files
    input: expand("results/filter/foreign/{{species}}.{fasta}.{{pool}}.eval", fasta=foreign_species)
    output: "results/filter/foreign/all.done"
    log: "logs/filter/foreign/all.log"
    conda: "../envs/base.yaml"
    shell: "touch {output} &> {log}"

max_mismatches = config["alignment_settings"]["max_mismatches"]
rule align_primers_to_foreign_species:
    input:
        expand("{index}/{{fasta}}.{version}.ebwt", version=range(1, 5), index=config["index_dir"]),
        expand("{index}/{{fasta}}.rev.{version}.ebwt", version=range(1, 3), index=config["index_dir"]),
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

rule format_align_primers_to_foreign_species:
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

rule dummy:
    input:
        "results/filter/foreign/{species}.{fasta}.{pool}.alignment.tsv",
        db = "results/{species}.db"
    output: "results/{species}.{fasta}.{pool}.eval"
    conda: "../envs/primers.yaml"
    shell: "touch {output}"
    