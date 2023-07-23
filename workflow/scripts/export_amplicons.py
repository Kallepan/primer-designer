"""
Extracts the discarded amplicons from the database and writes them to a json.
"""
import argparse
import json
import logging
import sys

import pandas as pd

from collections import defaultdict
from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extracts the discarded amplicons from the database and writes them to a json."
    )

    parser.add_argument(
        "--db",
        type=str,
        required=True,
        help="Path to the database",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )

    return parser.parse_args()


def __write_amplicons(amplicons: list[dict], output: str) -> None:
    """Write the amplicons to a json file."""
    with open(output, "w") as f:
        json.dump(amplicons, f, indent=4)


def __extract_failed_amplicons(db: DBHandler) -> pd.DataFrame:
    """Extract all failed amplicons from the database."""
    query = """
        SELECT amplicon_name, region_name, pool
        FROM failed_amplicons
    """

    data, columns = db.select(query)

    return pd.DataFrame(data, columns=columns)


def __extract_primers(db: DBHandler) -> pd.DataFrame:
    """Extract all proto_primers from the database."""
    query = """
        SELECT *
        FROM proto_primers
    """

    data, columns = db.select(query)

    return pd.DataFrame(data, columns=columns)


def __extract_amplicons_from_generation(df: pd.DataFrame) -> list[dict]:
    """
    Take the dataframe with all failed amplicons and iterate over them. Data is stored in a list with each item being an amplicon.
    """

    # Create a list of amplicons
    amplicons = []

    # Iterate over each amplicon
    for _, row in df.iterrows():
        amplicon = defaultdict(str)

        # Add the amplicon information
        amplicon["name"] = row["amplicon_name"]
        amplicon["region"] = row["region_name"]
        amplicon["pool"] = str(row["pool"])
        amplicon["discarded"] = True

        # Add the amplicon to the list
        amplicons.append(amplicon)

    return amplicons


def __extract_amplicons_from_filtering(df: pd.DataFrame) -> list[dict]:
    """
    Take the dataframe with all proto primers and group them by region. Data is stored in a list with each item being an amplicon.

    active_amplicons: list of amplicons that are active in the region
    discarded_amplicons: list of amplicons that are discarded in the region

    The amplicon is assigned each list by retrieving all prot_primers for an amplicon using groupby. If all forward or all reverse primers are discarded, the amplicon is discarded. Otherwise, it is active.
    """
    # Convert discarded column to boolean
    df["discarded"] = df["discarded"].astype(bool)

    # Group by region and amplicon
    dfs = df.groupby(["region_name", "amplicon_name", "pool"])

    # Create a list of amplicons
    amplicons = []

    # Iterate over each region and amplicon
    for (region_name, amplicon_name, pool), group in dfs:
        # group is a dataframe with all primers for an amplicon, region and pool
        # Split the dataframe into forward and reverse primers
        forward_primers = group[group["strand"] == "forward"]
        reverse_primers = group[group["strand"] == "reverse"]

        # Create a dictionary for the amplicon
        amplicon = defaultdict(str)

        # Add the amplicon information
        amplicon["name"] = amplicon_name
        amplicon["region"] = region_name
        amplicon["pool"] = str(pool)
        amplicon["n_forward"] = int(len(forward_primers))
        amplicon["n_reverse"] = int(len(reverse_primers))
        amplicon["n_discarded_forward"] = int(forward_primers["discarded"].sum())
        amplicon["n_discarded_reverse"] = int(reverse_primers["discarded"].sum())

        # If all forward or all reverse primers are discarded, the amplicon is discarded
        if forward_primers["discarded"].all() or reverse_primers["discarded"].all():
            amplicon["discarded"] = True
        else:
            amplicon["discarded"] = False

        # Add the amplicon to the list
        amplicons.append(amplicon)

    # Return the list of amplicons
    return amplicons


def main():
    args = __get_args()

    logging.info("Connecting to database...")
    db = DBHandler(args.db)

    # Extract the discarded amplicons generated during primer filtering
    logging.info("Extracting amplicons were primers were removed during filtering...")
    proto_primers = __extract_primers(db)
    amplicons = __extract_amplicons_from_filtering(proto_primers)

    # Extract the failed amplicons generated during primer generation
    logging.info("Extracting amplicons without primers...")
    failed_amplicons = __extract_failed_amplicons(db)
    amplicons.extend(__extract_amplicons_from_generation(failed_amplicons))

    # Write the amplicons to a json file
    __write_amplicons(amplicons, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
