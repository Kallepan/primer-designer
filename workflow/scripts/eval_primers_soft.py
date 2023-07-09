"""
Take each problematic primer (primer not aligning to the original site) and calculate a badness score for it.
This score consists of the following components:
    - Sum of 'Score' of each alignment (mismatches, mismatches_descriptor, matches) is precalculated in the database
    - If other complementary primers are found within a certain distance, the score is increased further
    - Optional hard filter to remove problematic primers with a score above a certain threshold
    - Soft filter (default) to calculate the score and store it in the database, but not remove the primers from the database
"""

import argparse
import logging
import sys

import pandas as pd

from handlers import DBHandler

DEFAULT_ADJACENCY_LIMIT = 500

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate the badness of the primers against target species"
    )
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
        "--species",
        type=str,
        required=True,
        help="Species to evaluate",
    )
    parser.add_argument(
        "--adjacency_limit",
        type=int,
        default=DEFAULT_ADJACENCY_LIMIT,
        help=f"Limit within which primers of opposite strand are considered adjancent enough to be problematic. Default: {DEFAULT_ADJACENCY_LIMIT}",
    )

    return parser.parse_args()


def __get_alignments_with_adjacent_primers(
    db: DBHandler, args: argparse.Namespace
) -> pd.DataFrame:
    """
    Returns all misaligned alignments with adjacent aligning primers
    I do a simple inner join with conditions to find adjacent alignments
    """
    data, columns = db.select(
        """
        WITH formatted_alignments AS (
            -- Select all alignments from the pool and add information from proto_primers
            SELECT 
                alignments.id, 
                alignments.primer_id, 
                alignments.position,
                alignments.matches,
                alignments.mismatches_descriptor,
                alignments.score, 
                alignments.aligned_to,
                alignments.pool,
                alignments.species,
                proto_primers.amplicon_name,
                proto_primers.strand AS primer_strand
            FROM alignments
            LEFT JOIN proto_primers
            ON proto_primers.id = alignments.primer_id
            WHERE 
                alignments.pool = ? AND 
                alignments.species = ?
        )
        SELECT
            alignments.id AS id,
            alignments.primer_id AS primer_id,
            alignments.position AS position,
            alignments.matches AS matches,
            alignments.mismatches_descriptor AS mismatches_descriptor,
            alignments.score AS score,
            alignments.aligned_to AS aligned_to,
            alignments.amplicon_name AS amplicon_name,
            alignments.pool AS pool,
            alignments.species AS species,
            adjacent_alignments.id AS adjacent_alignment_id,
            adjacent_alignments.primer_id AS adjacent_alignment_primer_id,
            adjacent_alignments.position AS adjacent_alignment_position,
            adjacent_alignments.score AS adjacent_alignment_score,
            adjacent_alignments.aligned_to AS adjacent_alignment_aligned_to,
            adjacent_alignments.amplicon_name AS adjacent_alignment_amplicon_name
        FROM formatted_alignments AS alignments
        -- inner join to only get alignments that have adjacent alignments
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
        ORDER BY 
            alignments.primer_id ASC, 
            alignments.id ASC
    """,
        (args.pool, args.species, args.adjacency_limit, args.adjacency_limit),
    )
    return pd.DataFrame(data, columns=columns)


def __calculate_badness_for_proto_primers(
    alignments: pd.DataFrame,
    adjacent_alignments: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate the badness score for each primer taking into account the alignments"""
    # sum all alignment scores for each primer
    primers_alignment_scores = (
        alignments.groupby("primer_id")["score"]
        .sum()
        .reset_index(name="alignment_score")
    )

    # sum all alignment scores for each adjacent primer for a primer
    adjacent_alignments_sum = (
        adjacent_alignments.groupby("primer_id")
        .agg(
            {
                "primer_id": "first",  # keep primer_id
                "score": "sum",
                "adjacent_alignment_score": "sum",
            }
        )
        .sum(axis=1)
        .reset_index(name="adjacency_score")
    )

    # merge dataframes on primer_id to prepare calculation of final badness score
    scores = pd.merge(
        primers_alignment_scores, adjacent_alignments_sum, on="primer_id", how="outer"
    )

    # Calculate badness score
    # TODO Implement calculation logic
    scores["badness"] = scores["alignment_score"] + scores["adjacency_score"]
    scores["badness"] = scores["badness"].fillna(0.0).astype(float)

    return scores


def __update_db_table(db: DBHandler, df: pd.DataFrame) -> None:
    """Updates the specified table with the provided dataframe"""
    db.executemany(
        f"""
        UPDATE proto_primers
        SET badness = ?
        WHERE id = ?
    """,
        df[["badness", "primer_id"]].values.tolist(),
    )


def __get_alignments(db: DBHandler, args: argparse.Namespace) -> pd.DataFrame:
    """Returns the specified table as a pandas dataframe"""
    data, columns = db.select(
        """
        SELECT * FROM alignments WHERE pool = ? AND species = ?
    """,
        (args.pool, args.species),
    )

    return pd.DataFrame(data, columns=columns)


def main():
    logging.info("Scoring primers...")
    args = __get_args()
    db = DBHandler(args.db)
    alignments = __get_alignments(db, args)
    alignments_with_adjacent_alignments = __get_alignments_with_adjacent_primers(
        db, args
    )
    scores = __calculate_badness_for_proto_primers(
        alignments, alignments_with_adjacent_alignments
    )

    __update_db_table(db, scores)

    # generate output file
    with open(args.output, "w") as f:
        scores.to_csv(f, sep="\t", index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
