import argparse
import logging
import pandas as pd

from db import DBHandler

DEFAULT_ADJACENCY_LIMIT = 500
DEFAULT_ALIGNMENTS_LIMIT = 3


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate the badness of the primers against target species. Hard Filter version"
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
    Returns all misaligned alignments with adjacent aligning primers for the target species (e.g.: missing correct alignments)
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
        -- Inner join with itself to get all adjacent alignments
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


def __get_primers_to_discard(
    args: argparse.Namespace, adjacent_alignments: pd.DataFrame
) -> tuple[list[str], pd.Series]:
    """
    marks the primers depending on their adjacent misaligned primers:
    Input: dataframe where each alignment is associated with its adjacent alignments
    Output: list of primers that should be marked as discarded

    Process:
        - Filter out the primers whose alignments do not meet the hard filter criteria
        - count the primers who still have mismatched alignments
        - Gather the ones above a threshold X
        - Mark them as discarded in the database
    """
    primer_adjacency_counts = adjacent_alignments["primer_id"].value_counts()
    primers_to_discard = primer_adjacency_counts[
        primer_adjacency_counts > DEFAULT_ALIGNMENTS_LIMIT
    ].index.tolist()  # TODO: add parameters
    return primers_to_discard, primer_adjacency_counts


def __update_db_tabke(db: DBHandler, primers_to_discard: list[str]) -> None:
    """Takes a list of primer_ids and marks them as discarded in the database"""
    db.executemany(
        f"""
        UPDATE proto_primers
        SET discarded = ?
        WHERE id = ?
    """,
        [(1, primer_id) for primer_id in primers_to_discard],
    )


def main():
    logging.info("Removing misaligned primers...")
    args = __get_args()
    db = DBHandler(args.db)
    adjacent_alignments = __get_alignments_with_adjacent_primers(db, args)
    primers_to_discard, counts = __get_primers_to_discard(args, adjacent_alignments)
    __update_db_tabke(db, primers_to_discard)

    # generate output file
    with open(args.output, "w") as f:
        counts.to_csv(f, sep="\t", index=True, header=True)

if __name__ == "__main__":
    main()
