import argparse
import subprocess
import sys

from io import StringIO
import pandas as pd
import sqlite3

DEFAULT_NUMBER_OF_MISMATCHES = 3

def get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align primers to reference genome and filter primers with multiple matches"
    )

    parser.add_argument(
        "--primers", type=str, required=True, help="Fasta file containing primers"
    )
    parser.add_argument(
        "--index", type=str, required=True, help="Path to index files"
    )
    parser.add_argument(
        "--db", type=str, required=False, help="Path to database file"
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=DEFAULT_NUMBER_OF_MISMATCHES,
        help=f"Number of mismatches allowed. Default: {DEFAULT_NUMBER_OF_MISMATCHES}",
    )

    return parser.parse_args()


def __run_bowtie(args: argparse.Namespace):
    shell_cmd = f"bowtie -v {args.mismatches} -a -x {args.index} -f {args.primers}"
    sp = subprocess.run(
        shell_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
    )
    if sp.returncode != 0:
        raise Exception(sp.stderr.decode("utf-8"))

    return sp.stdout.decode("utf-8")


def __parse_alignment(raw_alignment: str) -> pd.DataFrame:
    def decode_strand(symbol: str) -> str:
        if symbol == "+":
            return "forward"
        if symbol == "-":
            return "reverse"
        return "unknown"

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

    # Convert merged primer string to multiple columns
    # primer_region|primer_amplicon|primer_strand|primer_id
    # Reformat strand to forward/reverse
    # Drop primer column
    alignment[["primer_region", "primer_amplicon", "primer_strand", "primer_id"]] = alignment["primer"].apply(lambda x: pd.Series(str(x).split("|")))
    alignment["strand"] = alignment["strand"].apply(decode_strand)
    
    return alignment


def __write_to_csv(alignment: pd.DataFrame, output: str):
    alignment.to_csv(output, index=False, header=True, sep="\t")

def setup_db(args: argparse.Namespace) -> sqlite3.Connection:
    db = sqlite3.connect(args.db)

    # Create table if it doesn't exist
    # alignment consists of the following columns:
    db.execute(
        """
        CREATE TABLE IF NOT EXISTS primer_alignments ("""
    )


def main():
    print("Aligning primers to reference genome")
    args = get_parser()
    raw_alignment = __run_bowtie(args)
    alignment = __parse_alignment(raw_alignment)
    if args.output:
        __write_to_csv(alignment, args.output)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)