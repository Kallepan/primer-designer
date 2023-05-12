import argparse
import subprocess
import os
import sys

from io import StringIO
import pandas as pd


def get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align primers to reference genome and filter primers with multiple matches"
    )

    parser.add_argument(
        "--primers", type=str, required=True, help="Fasta file containing primers"
    )
    parser.add_argument("--index", type=str, required=True, help="Path to index files")
    parser.add_argument("--output", type=str, help="Output file name")

    return parser.parse_args()


def __run_bowtie(args: argparse.Namespace):
    shell_cmd = f"bowtie -v 0 -a -x {args.index} -f {args.primers}"
    sp = subprocess.run(
        shell_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
    )
    if sp.returncode != 0:
        raise Exception(sp.stderr.decode("utf-8"))

    return sp.stdout.decode("utf-8")


def __parse_alignment(raw_alignment: str):
    """Takes the tab seperated output from the alignment and parses it to a pandas dataframe"""
    tsv_string = StringIO(raw_alignment)

    alignment = pd.read_csv(
        tsv_string,
        sep="\t",
        header=None,
        names=[
            "primer",
            "strand",
            "chromosome",
            "position",
            "sequence",
            "read_quality",
            "matches",
            "mismatches_descriptor",
        ],
    )
    # matches is reported as additional matches, therefore we need to add 1 to get the actual number of matches
    alignment["matches"] = alignment["matches"].apply(lambda x: int(x) + 1)
    return alignment


def __write_to_csv(alignment: pd.DataFrame, output: str):
    alignment.to_csv(output, index=False, header=True, sep=",")


def main():
    args = get_parser()
    raw_alignment = __run_bowtie(args)
    alignment = __parse_alignment(raw_alignment)
    if args.output:
        __write_to_csv(alignment, args.output)


if __name__ == "__main__":
    main()
