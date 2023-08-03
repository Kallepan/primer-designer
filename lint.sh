#! /bin/bash

# the script is meant to lint all rules and scripts in a snakemake workflow before push

echo "##########################################################################################"
echo "execute black to lint python scripts"
echo "##########################################################################################"

# To find the location of the commands origin, type "which <command>" in the command-line and set
# the path in "PATH=$PATH:/bin:<result_path>"
# For example: "which black" -> /homes/kkandeepan/.conda/envs/sm/bin/black
# so the line here has to be PATH=$PATH:/bin:/homes/kkandeepan/.conda/envs/sm/bin/black
# and the Path needs to be exported

PATH=$PATH:/bin:/homes/kkandeepan/mambaforge/envs/snakemake/bin/black
export PATH

dir="workflow/scripts/*.py"
for i in $dir; do echo $i; black $i;done

dir="workflow/scripts/primer_gen/*.py"
for i in $dir; do echo $i; black $i;done

echo "##########################################################################################"
echo "execute snakemake linter to lint rules"
echo "##########################################################################################"

PATH=$PATH:/bin:/homes/kkandeepan/mambaforge/envs/snakemake/bin/snakemake
export PATH

snakemake --lint

