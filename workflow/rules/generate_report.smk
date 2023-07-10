rule build_report:
    conda: "../envs/nodejs.yaml"
    input: 
        primer_set = "results/final/{species}.primer_set.json",
        loss = "results/final/{species}.loss_set.json",
        regions = "results/{species}.regions.json"
    log: "logs/report/{species}.build_report.log"
    output: "results/final/{species}.summary.html"
    params:
        source_dir = "packages/visualizer/"
    shell:
        """
            bash workflow/scripts/build_report.sh {input.primer_set} {input.loss} {input.regions} {params.source_dir} {output} &> {log}
        """
