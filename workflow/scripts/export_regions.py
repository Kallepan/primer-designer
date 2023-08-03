import argparse
import logging

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

    parser.add_argument(
        "--plotting_buffer",
        type=int,
        default=300,
        required=False,
        help="Buffer to add to the start and end positions of the regions to be plotted",
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


def __extract_regions_sequence_from_fasta(
    args: argparse.Namespace, region_df: pd.DataFrame
) -> pd.DataFrame:
    """Extract the regions sequence from the fasta file and add them to the region dataframe"""

    def __extract_region_from_seqrecord(
        start: int, end: int, seqrecord: SeqIO.SeqRecord
    ) -> str:
        """Extract the region from the seqrecord"""
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


def __adjust_region_df_for_plotting(
    args: argparse.Namespace, region_df: pd.DataFrame
) -> pd.DataFrame:
    def __adjust_start(row: pd.Series) -> int:
        """Adjust the start position of the region to be plotted"""
        return max(row["start"] - args.plotting_buffer, 0)

    def __adjust_end(row: pd.Series, seqrecord: SeqIO.SeqRecord) -> int:
        """Adjust the end position of the region to be plotted"""
        return min(row["end"] + args.plotting_buffer, len(seqrecord))

    def __extract_region_from_seqrecord(
        start: int, end: int, seqrecord: SeqIO.SeqRecord
    ) -> str:
        """Extract the region from the seqrecord"""
        return str(seqrecord.seq[start:end])

    # Load in the fasta file as seqrecord
    with open(args.fasta, "r") as file:
        seqrecord = SeqIO.read(file, "fasta")

    # Adjust the start and end positions of the regions
    region_df["start"] = region_df.apply(__adjust_start, axis=1)
    region_df["end"] = region_df.apply(lambda row: __adjust_end(row, seqrecord), axis=1)

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
    region_df = __extract_regions_sequence_from_fasta(args, region_df)

    __to_db(region_df, db)

    logging.info("Writing regions to json file")
    region_df = __adjust_region_df_for_plotting(args, region_df)
    __to_json(region_df, args)


if __name__ == "__main__":
    main()
