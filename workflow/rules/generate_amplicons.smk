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


""" Split amplicons into separate files. """
rule split_by_amplicon:
    input:
        amplicons_fasta = "results/{sample}.amplicons.fasta"
    log:
        "logs/split_by_amplicon.{sample}.log"
    conda:
        "../envs/amplicons.yaml"
    output:
        "results/{sample}.split"
    shell:
        """cat {input.amplicons_fasta} | awk '{{
            if (substr($0, 1, 1)==">") {{ filename=("results/amplicons/{wildcards.sample}." substr($0, 2) ".fasta") }}
            print $0 >> filename
            close(filename)
            }}'

            if [ $? -eq 0 ]; then
                echo "done" > {output}
            else
                echo "Failed to split amplicons" > {log}
            fi
        """ 