import logging
import argparse
import pandas as pd

from db import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Formatting primer alignments")

    parser.add_argument("--input", type=str, required=True, help="Path to input file")
    parser.add_argument("--output", type=str, required=True, help="Path to output file")
    parser.add_argument("--pool", type=str, required=True, help="Pool number")
    parser.add_argument("--species", type=str, required=True, help="Species name")
    parser.add_argument("--db", type=str, required=True, help="Path to database file")

    return parser.parse_args()


def __parse_alignment(args: argparse.Namespace) -> pd.DataFrame:
    def decode_aligned_to(symbol: str) -> str:
        if symbol == "+":
            return "forward"
        if symbol == "-":
            return "reverse"
        return "unknown"

    """Takes the tab seperated output from the alignment and parses it to a pandas dataframe"""

    alignment = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        names=[
            "primer_id",
            "aligned_to",
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

    # Reformat aligned_to to forward/reverse
    alignment["aligned_to"] = alignment["aligned_to"].apply(decode_aligned_to)
    alignment["pool"] = args.pool
    alignment["species"] = args.species
    return alignment


def main():
    logging.info("Formatting primer alignments")
    args = __get_parser()
    db = DBHandler(args.db)

    # read in alignment
    alignment = __parse_alignment(args)

    # write output to csv and database
    alignment.to_sql("alignments", db.con, if_exists="append", index=False)

    alignment.to_csv(args.output, index=False)

    logging.info("Wrote primer alignments to database")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(e)
        raise e
