# Global configuration
pools = range(0, config["metadata"]["num_pools"])


rule all:
    input:
        expand("{results}/dump/{{species}}.proto_primers.tsv", results=config["results_dir"]),
        expand("{results}/dump/{{species}}.alignments.tsv", results=config["results_dir"]),
        expand("{results}/{{species}}.bed", results=config["results_dir"]),
        expand("{results}/{{species}}.summary.html", results=config["results_dir"]),
        expand("{results}/{{species}}.primer_set.tsv", results=config["results_dir"])
    output: touch("results/{species}.dummy")