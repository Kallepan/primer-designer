# Fetch all fasta files to be excluded against
import os
fasta_files = os.listdir("data/exclude/")

rule all_foreign_species:
    input:
        expand("results/filter/foreign/{{species}}.{fasta}.{{pool}}.eval", fasta=fasta_files)
    output: "results/filter/foreign/all.done"
    log: "logs/filter/foreign/all.log"
    conda: "../envs/base.yaml"
    shell: "touch {output} &> {log}"

# Build indexes. Note: Indexes can be supplied to skip this step
rule build_indexes_for_foreign_species:
    input: 
        foreign_species_fasta = "data/exclude/{fasta}.fasta"
    output:
        expand("data/indexes/{{fasta}}.{version}.ebwt", version=range(1, 5)),
        expand("data/indexes/{{fasta}}.rev.{version}.ebwt", version=range(1, 3)),
    log: "logs/indexes/{fasta}.log"
    params:
        outdir = lambda w, input: os.path.join("data", "indexes", os.path.splitext(os.path.basename(input.foreign_species_fasta))[0])
    conda: "../envs/bowtie.yaml"
    shell: "bowtie-build {input.foreign_species_fasta} {params.outdir} &> {log}"

max_mismatches = config["alignment_settings"]["max_mismatches"]
rule align_primers_to_foreign_species:
    input:
        expand("data/indexes/{{fasta}}.{version}.ebwt", version=range(1, 5)),
        expand("data/indexes/{{fasta}}.rev.{version}.ebwt", version=range(1, 3)),
        primers_fasta="results/{species}.{pool}.proto_primers.fasta",
        db = "results/{species}.db"
    output: "results/filter/foreign/{species}.{fasta}.{pool}.alignment.raw"
    log: "logs/filter/foreign/{species}.{fasta}.{pool}.alignment.log"
    params:
        index = lambda w: os.path.join("data", "indexes", w.fasta),
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
    