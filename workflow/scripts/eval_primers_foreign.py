import logging
import argparse
import sys
import pandas as pd

from handlers import DBHandler

logging.basicConfig(level=logging.INFO)

DEFAULT_ADJACENCY_LIMIT = 500


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate alignment of primers to foreign genomes"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    parser.add_argument(
        "--db",
        type=str,
        required=True,
        help="Path to the database file",
    )
    parser.add_argument(
        "--pool",
        type=str,
        required=True,
        help="Pool number to evaluate",
    )
    parser.add_argument(
        "--adjacency_limit",
        type=int,
        default=DEFAULT_ADJACENCY_LIMIT,
        help=f"Limit within which primers of opposite strand are considered adjancent enough to be problematic. Default: {DEFAULT_ADJACENCY_LIMIT}",
    )

    return parser.parse_args()


def main():
    logging.info("Scoring primers against foreign genomes")
    args = __get_args()

    db = DBHandler(args.db)

    # Get alignments score and adjacent alignments score
    # alignments =


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
