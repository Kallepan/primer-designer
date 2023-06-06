"""
Take each problematic primer (primer not aligning to the original site) and calculate a badness score for it.
This score consists of the following components:
    - Sum of 'Score' of each alignment (mismatches, mismatches_descriptor, matches) is precalculated in the database
    - If other complementary primers are found within a certain distance, the score is increased further
    - Optional hard filter to remove problematic primers with a score above a certain threshold
    - Soft filter (default) to calculate the score and store it in the database, but not remove the primers from the database
"""

import argparse
import pandas as pd

from db import DBHandler

DEFAULT_ADJACENCY_LIMIT = 500

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate the badness of the primers against target species")
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="JSON File with primer sequences and badness score",
    )
    parser.add_argument(
        "--db",
        type=str,
        required=True,
        help="Path to database file containing primer sequences and alignments",
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
    parser.add_argument(
    # TODO: implement hard filter
        "--hard_filter",
        action="store_true",
        help="If set, primers with a badness score above a certain threshold are marked as discarded in the database",
    )

    return parser.parse_args()

def get_alignments_with_adjacent_primers(db: DBHandler, args: argparse.Namespace) -> pd.DataFrame:
    """ Returns all alignments with adjacent primers """
    data, columns = db.select("""
        WITH formatted_alignments AS (
            -- Select all alignments from the pool and add information from proto_primers
            SELECT 
                alignments.id, 
                alignments.primer_id, 
                alignments.position, 
                alignments.score, 
                alignments.aligned_to,
                alignments.pool,
                proto_primers.amplicon_name
            FROM alignments
            LEFT JOIN proto_primers
            ON 
                proto_primers.id = alignments.primer_id AND
                proto_primers.pool = alignments.pool
            WHERE alignments.pool = ?
        )
        SELECT
            alignments.id AS alignment_id,
            alignments.primer_id AS alignment_primer_id,
            alignments.position AS alignment_position,
            alignments.score AS alignment_score,
            alignments.aligned_to AS alignment_aligned_to,
            alignments.amplicon_name AS alignment_amplicon_name,
            adjacent_alignments.id AS adjacent_alignment_id,
            adjacent_alignments.primer_id AS adjacent_alignment_primer_id,
            adjacent_alignments.position AS adjacent_alignment_position,
            adjacent_alignments.score AS adjacent_alignment_score,
            adjacent_alignments.aligned_to AS adjacent_alignment_aligned_to,
            adjacent_alignments.amplicon_name AS adjacent_alignment_amplicon_name,
            alignments.pool AS pool
        FROM formatted_alignments AS alignments
        INNER JOIN formatted_alignments AS adjacent_alignments
        ON 
            -- Select all alignments from the pool that are adjacent to the current alignment. Ignore the same amplicon and the same strand
            adjacent_alignments.amplicon_name <> alignments.amplicon_name AND 
            adjacent_alignments.aligned_to <> alignments.aligned_to AND
            CASE
                -- Find only adjacent alignments (withing adjacency_limit) on the same strand
                WHEN alignments.aligned_to = 'forward' THEN
                    adjacent_alignments.position >= alignments.position AND
                    adjacent_alignments.position <= alignments.position + ?
                WHEN alignments.aligned_to = 'reverse' THEN
                    adjacent_alignments.position <= alignments.position AND
                    adjacent_alignments.position >= alignments.position - ?
            END
        ORDER BY alignments.id ASC
    """, (args.pool, args.adjacency_limit, args.adjacency_limit))
    return pd.DataFrame(
        data,
        columns=columns
    )

def calculate_badness_for_proto_primers(args: argparse.Namespace, proto_primers: pd.DataFrame, alignments: pd.DataFrame, adjacent_alignments: pd.DataFrame) -> pd.DataFrame:
    """ Calculate the badness score for each primer taking into account the alignments """
    # sum all alignment scores for each primer
    primers_alignment_scores = alignments.groupby("primer_id")["score"].sum()
    
    # sum all alignment scores for each primer with adjacent primers

def filter_primers_with_multiple_alignments(alignments: pd.DataFrame) -> pd.DataFrame:
    pass

def get_db_table(db: DBHandler, table_name: str, args: argparse.Namespace) -> pd.DataFrame:
    """ Returns the specified table as a pandas dataframe """
    data, columns = db.select(f"SELECT * FROM {table_name} WHERE pool = ?", (args.pool,))
    return pd.DataFrame(
        data,
        columns=columns
    )

def main():
    print("Evaluating primers against target species")
    args = get_args()
    db = DBHandler(args.db)
    proto_primers = get_db_table(db, "proto_primers", args)
    alignments = get_db_table(db, "alignments", args)
    alignments_with_adjacent_alignments = get_alignments_with_adjacent_primers(db, args)
    alignments_with_adjacent_alignments.to_csv("temp.csv", index=False)
    calculate_badness_for_proto_primers(args, proto_primers, alignments, alignments_with_adjacent_alignments)

if __name__ == "__main__":
    main()