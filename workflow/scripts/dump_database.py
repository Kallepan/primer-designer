"""
Script to dump the database to a set of TSV files for each table in the database.
This is useful for debugging and for exporting the database to a different format.
"""
import argparse
import logging
import pandas as pd
import os

from handlers import DBHandler
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dump the results from the database")
    parser.add_argument("--db", type=str, help="Path to the database", required=True)
    parser.add_argument("--species", type=str, help="Species name", required=True)
    parser.add_argument(
        "--output_dir", type=str, help="Path to the alignment output dir", required=True
    )

    return parser.parse_args()


def __export(args: argparse.Namespace) -> None:
    db = DBHandler(args.db)
    tables = db.get_tables()

    # create output dir if it does not exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    for table_name in tables:
        # Extract the tables
        path = os.path.join(args.output_dir, f"{args.species}.{table_name}.tsv")
        table_data, columns = db.select(f"SELECT * FROM {table_name}")
        df = pd.DataFrame(table_data, columns=columns)
        df.to_csv(path, sep="\t", index=False)


def main():
    logging.info("Dumping the database")
    args = __get_args()
    __export(args)


if __name__ == "__main__":
    main()
