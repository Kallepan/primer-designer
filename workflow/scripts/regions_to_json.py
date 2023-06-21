import argparse
import logging
import sys

from db import DBHandler
import pandas as pd

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format the regions into a json file to be used during report generation"
    )
    parser.add_argument(
        "--regions", type=str, required=True, help="Path to the regions file"
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

    return df


def __to_db(regions: pd.DataFrame, db: DBHandler) -> None:
    regions.to_sql("regions", db.con, if_exists="append", index=False)


def __to_json(region_df: pd.DataFrame, args: argparse.Namespace) -> None:
    region_df.to_json(args.output, orient="records", indent=4)


def main():
    logging.info("Formatting regions into a json file and writing to database")
    args = __get_args()
    db = DBHandler(args.db)
    region_df = __get_regions(args.regions)
    __to_db(region_df, db)
    __to_json(region_df, args)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
