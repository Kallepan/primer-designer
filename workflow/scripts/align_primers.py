import argparse
import subprocess
import sys

from io import StringIO
from db import DBHandler
import pandas as pd

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
        "--output", type=str, required=False, help="Path to output file"
    )
    parser.add_argument(
        "--pool", type=str, required=False, help="Pool number"
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


def __parse_alignment(raw_alignment: str, args: argparse.Namespace) -> pd.DataFrame:
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
            "primer_id",
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

    # Reformat strand to forward/reverse
    alignment["strand"] = alignment["strand"].apply(decode_strand)
    alignment["pool"] = args.pool
    return alignment

def main():
    print("Aligning primers to reference genome")

    args = get_parser()
    db = DBHandler(args.db)
    db.setup_alignments_table()
    raw_alignment = __run_bowtie(args)
    alignment = __parse_alignment(raw_alignment, args)

    # write output to csv and database
    alignment.to_sql("alignments", db.con, if_exists="append", index=False, chunksize=1000)
    alignment.to_csv(args.output, index=False)
    
    print("Wrote primer alignments to database")

if __name__ == "__main__":
    main()
    try:
        pass
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)