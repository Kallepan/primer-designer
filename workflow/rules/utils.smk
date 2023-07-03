rule create_db:
    input: 
        regions_sql_file = "workflow/sql/regions.sql",
        proto_primers_sql_file = "workflow/sql/proto_primers.sql",
        alignments_sql_file = "workflow/sql/alignments.sql",
    output: "results/{species}.db"
    log: "logs/{species}.db.log"
    conda: "../envs/dump.yaml"
    shell:
        """
        python3 workflow/scripts/create_db.py \
            --db {output} &> {log} \
            --regions_sql_file {input.regions_sql_file} \
            --proto_primers_sql_file {input.proto_primers_sql_file} \
            --alignments_sql_file {input.alignments_sql_file} \
            &> {log}
        """
