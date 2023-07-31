rule build_report:
    conda: "../envs/nodejs.yaml"
    input: 
        primer_set = expand("{results}/{{species}}.primer_set.json", results = config["results_dir"]),
        loss = expand("{results}/{{species}}.loss_set.json", results = config["results_dir"]),
        regions = expand("{results}/{{species}}.regions.json", results = config["results_dir"]),
        amplicons = expand("{results}/{{species}}.amplicons.json", results = config["results_dir"])
    log: "logs/report/{species}.build_report.log"
    output: expand("{results}/{{species}}.summary.html", results = config["results_dir"])
    shell:
        """
        bash workflow/scripts/build_report.sh \
            {input.primer_set} \
            {input.loss} \
            {input.regions} \
            {input.amplicons} \
            {output} \
            &> {log}
        """

