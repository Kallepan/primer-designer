"""
Creates the sqlite database with the tables:
    - regions
    - alignments
    - proto_primers
See the sql files in the sql directory for more information.
"""
import argparse
import logging

from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create the sqlite database")

    parser.add_argument(
        "--db", type=str, required=True, help="Path to the sqlite database"
    )
    parser.add_argument(
        "--regions_sql_file",
        type=str,
        required=True,
        help="Path to the regions.sql file",
    )
    parser.add_argument(
        "--alignments_sql_file",
        type=str,
        required=True,
        help="Path to the alignments.sql file",
    )
    parser.add_argument(
        "--proto_primers_sql_file",
        type=str,
        required=True,
        help="Path to the proto_primers.sql file",
    )

    return parser.parse_args()


def main():
    logging.info("Creating the sqlite database")

    args = __get_args()
    db = DBHandler(path_to_db=args.db)
    db.setup_table(args.regions_sql_file)
    db.setup_table(args.proto_primers_sql_file)
    db.setup_table(args.alignments_sql_file)


if __name__ == "__main__":
    main()
