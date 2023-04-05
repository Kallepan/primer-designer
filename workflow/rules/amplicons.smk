rule:
    input:
        fasta_file = "data/{sample}.fasta",
        loci_file = "data/loci_formatted.csv"
    output:
        output_file = "results/{sample}.amplicons.fasta"
    conda:
        "../envs/amplicons.yaml"
    params:
        amplicon_size = config["amplicon_size"],
        amplicon_buffer = config["amplicon_buffer_size"]
    log:
        out = "logs/amplicons_{sample}.out",
        err = "logs/amplicons_{sample}.err"
    shell:
        "python workflow/scripts/generate_amplicons.py {input.fasta_file} {input.loci_file} {output.output_file} -s {params.amplicon_size} -b {params.amplicon_buffer} 2> {log.err} > {log.out}"