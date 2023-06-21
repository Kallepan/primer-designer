"""
    Use alignment with target species to determin the position of the primers in the region
"""
import logging
import argparse
import sys

import pandas as pd
from db import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Determine the position of the primers in the region"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="TSV dump of the primers with the position in the region",
    )

    parser.add_argument(
        "--db",
        type=str,
        required=True,
        help="Path to database file containing primer sequences and alignments",
    )

    # Necessary to get the correct alignments from the databse
    parser.add_argument(
        "--species",
        type=str,
        required=True,
        help="Species to evaluate",
    )
    parser.add_argument(
        "--pool",
        type=str,
        required=True,
        help="Pool number to evaluate",
    )

    return parser.parse_args()


def __get_correct_alignments_with_primers(
    args: argparse.Namespace, db: DBHandler
) -> pd.DataFrame:
    alignments, columns = db.select(
        """
        SELECT proto_primers.id, alignments.position
        FROM proto_primers
        LEFT JOIN alignments
        ON proto_primers.id = alignments.primer_id
        WHERE
            alignments.pool = ? AND
            alignments.species = ? AND (
                -- filter setting to get the correct alignments, e.g.: ignore misalignments
                alignments.matches = 1 AND
                alignments.mismatches_descriptor IS NULL AND
                proto_primers.strand = alignments.aligned_to 
            )
        ORDER BY alignments.id ASC
        
    """,
        (args.pool, args.species),
    )

    return pd.DataFrame(alignments, columns=columns, dtype="Int64")


def __update_db_with_positions(
    args: argparse.Namespace, primer_df: pd.DataFrame, db: DBHandler
) -> None:
    # First check if each primer has only one alignment
    if primer_df["id"].nunique() != primer_df.shape[0]:
        primer_df.to_csv(args.output, sep="\t", index=False)
        exception = Exception("Something went wrong: Some primers have more than one alignment. Check the output file.")
        logging.exception(exception)
        raise exception

    # Update the database with the position of the primers
    db.executemany(
        """
        UPDATE proto_primers
        SET position = ?
        WHERE id = ?
    """,
        primer_df[["position", "id"]].values.tolist(),
    )

    # Write positions to file
    primer_df.to_csv(args.output, sep="\t", index=False)


def main():
    logging.info("Calculating primer positions")

    args = __get_args()
    db = DBHandler(args.db)
    primer_df = __get_correct_alignments_with_primers(args, db)
    __update_db_with_positions(args, primer_df, db)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
