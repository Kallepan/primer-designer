import argparse
import logging
import sys

from db import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create the sqlite database")

    parser.add_argument(
        "--db", type=str, required=True, help="Path to the sqlite database"
    )

    return parser.parse_args()


def main():
    logging.info("Creating the sqlite database")

    args = __get_args()
    db = DBHandler(path_to_db=args.db)
    db.setup_proto_primers_table()
    db.setup_alignments_table()
    db.setup_regions_table()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
