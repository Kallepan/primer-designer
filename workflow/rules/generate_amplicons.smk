rule generate_amplicons:
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
        out = "logs/{sample}.amplicons.log"
    shell:
        """
        python workflow/scripts/generate_amplicons.py {input.fasta_file} {input.loci_file} {output.output_file} -s {params.amplicon_size} -b {params.amplicon_buffer} > {log.out};
        
        # separate amplicons into individual files
        cat {output.output_file} | awk '{{
            if (substr($0, 1, 1)==">") {{ filename=("results/amplicons/{wildcards.sample}." substr($0, 2) ".fasta") }}
            print $0 >> filename
            close(filename)
            }}'

            if [ $? -eq 0 ]; then
                echo "Amplicons split" > {log.out}
            else
                echo "Failed to split amplicons" > {log.out}
            fi
        """