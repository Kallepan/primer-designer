import argparse
import logging
import sys

from handlers import DBHandler
import pandas as pd
import numpy as np
from Bio import SeqIO

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format the regions into a json file to be used during report generation"
    )
    parser.add_argument(
        "--regions", type=str, required=True, help="Path to the regions file"
    )
    parser.add_argument(
        "--fasta", type=str, required=True, help="Path to the fasta file"
    )
    parser.add_argument(
        "--db", type=str, required=True, help="Path to the sqlite database"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to the output json file"
    )

    return parser.parse_args()


def __get_regions(path_to_file: str) -> pd.DataFrame:
    df = pd.read_csv(
        path_to_file,
        sep=",",
        header=0,
        dtype={
            "name": "string",
            "start": "int64",
            "end": "int64",
        },
    )

    # switch up the start and end columns if the start is greater than the end
    df["start"], df["end"] = np.where(
        df["start"] > df["end"], (df["end"], df["start"]), (df["start"], df["end"])
    )
    return df


def __extract_amplicons_from_fasta(
    args: argparse.Namespace, region_df: pd.DataFrame
) -> pd.DataFrame:
    """Extract the amplicons from the fasta file and add them to the region dataframe"""

    def __extract_region_from_seqrecord(
        start: int, end: int, seqrecord: SeqIO.SeqRecord
    ) -> str:
        return str(seqrecord.seq[start:end])

    # Load in the fasta file as seqrecord
    with open(args.fasta, "r") as file:
        seqrecord = SeqIO.read(file, "fasta")

    # Extract the regions from the fasta file
    region_df["sequence"] = region_df.apply(
        lambda row: __extract_region_from_seqrecord(
            row["start"], row["end"], seqrecord
        ),
        axis=1,
    )

    return region_df


def __to_db(regions: pd.DataFrame, db: DBHandler) -> None:
    # check if regions table has entries
    rows, _ = db.select("SELECT * FROM regions")
    if len(rows) > 0:
        logging.warning("Regions table already has entries, skipping")
        return

    # Append to regions table
    regions.to_sql("regions", db.conn, if_exists="append", index=False)


def __to_json(region_df: pd.DataFrame, args: argparse.Namespace) -> None:
    region_df.to_json(args.output, orient="records", indent=4)


def main():
    logging.info("Formatting regions into a json file and writing to database")
    args = __get_args()
    db = DBHandler(args.db)
    region_df = __get_regions(args.regions)
    region_df = __extract_amplicons_from_fasta(args, region_df)
    __to_db(region_df, db)
    __to_json(region_df, args)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
