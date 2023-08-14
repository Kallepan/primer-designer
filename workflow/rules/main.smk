# Global configuration
pools = range(0, config["metadata"]["pool_count"])


rule all:
    input:
        "results/{species}/dump/.done",
        expand("{results}/{{species}}.bed", results=config["results_dir"]),
        expand("{results}/{{species}}.summary.html", results=config["results_dir"]),
        expand("{results}/{{species}}.primer_set.tsv", results=config["results_dir"])
    output: touch("results/{species}.dummy")