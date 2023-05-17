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
        "--output", type=str, required=False, help="Path to output file"
    )
    parser.add_argument(
        "--pool", type=int, required=False, help="Pool number"
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

def setup_db(args: argparse.Namespace) -> sqlite3.Connection:
    db = sqlite3.connect(args.db)

    # Create table if it doesn't exist
    # alignment consists of the following columns:
    # pool|primer_id|strand|chromosome|position|sequence|read_quality|matches|mismatches_descriptor
    # primer_id is a foreign key to the proto_primers table
    db.execute(
        """
        CREATE TABLE IF NOT EXISTS primer_alignments (
            pool INT NOT NULL,
            primer_id TEXT NOT NULL,
            strand TEXT NOT NULL,
            chromosome TEXT NOT NULL,
            position INT NOT NULL,
            sequence TEXT NOT NULL,
            read_quality TEXT NOT NULL,
            matches INT NOT NULL,
            mismatches_descriptor TEXT,

            FOREIGN KEY (primer_id) REFERENCES proto_primers (primer_id)
        );
        """
    )

    db.execute(
        """
        CREATE INDEX IF NOT EXISTS primer_alignments_index ON primer_alignments (pool, primer_id, position, matches);
        """
    )

    db.commit()
    return db

def main():
    print("Aligning primers to reference genome")

    args = get_parser()
    db = setup_db(args)
    raw_alignment = __run_bowtie(args)
    alignment = __parse_alignment(raw_alignment, args)

    # write output to csv and database
    alignment.to_sql("primer_alignments", db, if_exists="append", index=False)
    alignment.to_csv(args.output, index=False)
    
    db.commit()
    print("Wrote primer alignments to database")
    
    db.close()

if __name__ == "__main__":
    main()
    try:
        pass
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)