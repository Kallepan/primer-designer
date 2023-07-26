import argparse
import logging
import sys

import pandas as pd

from collections import defaultdict
from handlers import DBHandler, Graph

logging.basicConfig(level=logging.INFO)

DEFAULT_ADJACENCY_LIMIT = 500
DEFAULT_BASES_TO_IGNORE = 3
DEFAULT_MAX_MISALIGNMENTS = -1


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluates the primers against the original target species"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output file",
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
        help="Pool id to evaluate",
    )
    parser.add_argument(
        "--species",
        type=str,
        required=True,
        help="Species to evaluate",
    )

    # Optional Args
    parser.add_argument(
        "--adjacency_limit",
        type=int,
        default=DEFAULT_ADJACENCY_LIMIT,
        help=f"Limit within which primers of opposite strand are considered adjancent enough to be problematic. Default: {DEFAULT_ADJACENCY_LIMIT}",
    )
    parser.add_argument(
        "--bases_to_ignore",
        type=int,
        default=DEFAULT_BASES_TO_IGNORE,
        help=f"Number of bases with mismatches at the 3' end of the primer to ignore when evaluating alignments. Default: {DEFAULT_BASES_TO_IGNORE}",
    )
    parser.add_argument(
        "--max_misalignments",
        type=int,
        default=DEFAULT_MAX_MISALIGNMENTS,
        help=f"Number of misalignments allowed in the primer alignment. Default: -1. If set to {DEFAULT_MAX_MISALIGNMENTS}, no limit is applied.",
    )

    return parser.parse_args()


def __get_alignments(db: DBHandler, args: argparse.Namespace) -> pd.DataFrame:
    """Fetches all alignments from a given species and pool."""
    # Build the query
    query = """
    SELECT
        alignments.id,
        alignments.primer_id, 
        alignments.position,
        alignments.matches, 
        alignments.mismatches_descriptor,
        alignments.aligned_to, 
        alignments.chromosome, 
        proto_primers.amplicon_name, 
        proto_primers.strand AS primer_strand,
        LENGTH(proto_primers.sequence) AS primer_length
    FROM alignments
    LEFT JOIN proto_primers ON proto_primers.id = alignments.primer_id
    WHERE 
        alignments.species = ? AND 
        alignments.pool = ? AND 
        proto_primers.discarded = 0 AND (
        --filter out all alignments that are correct
        alignments.matches <> 1 OR
        alignments.mismatches_descriptor IS NOT NULL OR
        proto_primers.strand <> alignments.aligned_to OR
        proto_primers.position <> alignments.position
    ) AND ((
        mismatches_descriptor IS NOT NULL OR
        LENGTH(proto_primers.sequence) < 10
    ) OR (
    """

    # Add the bases to ignore
    for i in range(1, args.bases_to_ignore + 1):
        query += f""" INSTR(alignments.mismatches_descriptor, CAST(LENGTH(proto_primers.sequence) - {i} AS TEXT)) = 0 AND """

    # Finish the query
    query += """
        1 = 1 ))
    ORDER BY alignments.chromosome ASC, alignments.position ASC
    """

    alignments, columns = db.select(query, (args.species, args.pool))

    return pd.DataFrame(alignments, columns=columns)


def __calculate_adjacent_alignments(
    alignments: pd.DataFrame,
    adjacency_limit: int,
) -> dict[str, list]:
    """
    For each alignment, find the next alignment that is within the adjacency limit.
    Keep track of the adjacent alignments for each primer in a dict.
    """

    # Split alignments into tables for each chromosome
    alignments_by_chromosome = alignments.groupby("chromosome")

    adjacent_primers_for_primer: dict[str, list] = defaultdict(list)
    for chromosome, alignments in alignments_by_chromosome:
        logging.info(
            f"Processing alignments for chromosome {chromosome} with {len(alignments)} alignments"
        )

        # Sort alignments by position
        alignments = alignments.sort_values("position")

        # Split the alignments into forward and reverse alignments
        alignments_forward = alignments[alignments["aligned_to"] == "forward"]
        alignments_reverse = alignments[alignments["aligned_to"] == "reverse"]

        # For each alignment, in teh dataframes, find all alignments of the complementary strand that are within the adjacency limit
        for _, alignment in alignments_forward.iterrows():
            adjacent_alignments = alignments_reverse[
                (alignments_reverse["position"] >= alignment["position"])
                & (
                    alignments_reverse["position"]
                    <= alignment["position"] + adjacency_limit
                )
            ]

            # if adjacent alignments were empty continue
            if adjacent_alignments.empty:
                continue

            # Turn dataframe into a list of dicts
            adjacent_alignments_info = adjacent_alignments.to_dict("records")

            # Add the adjacent alignments to the dict
            key = int(alignment["primer_id"])
            adjacent_primers_for_primer[key].extend(adjacent_alignments_info)

        # Repeat for reverse alignments
        for _, alignment in alignments_reverse.iterrows():
            # Filter out complementary alignments that are not within the adjacency limit
            adjacent_alignments = alignments_forward[
                (
                    alignments_forward["position"]
                    >= alignment["position"] - adjacency_limit
                )
                & (alignments_forward["position"] <= alignment["position"])
            ]

            # if adjacent_alignments is empty, we can continue
            if adjacent_alignments.shape[0] == 0:
                continue

            # Turn dataframe into a list of dicts
            adjacent_alignments_info = adjacent_alignments.to_dict("records")

            # Add the adjacent alignments to the dict
            key = int(alignment["primer_id"])
            adjacent_primers_for_primer[key].extend(adjacent_alignments_info)

    return adjacent_primers_for_primer


def __select_misaligned_primers(
    alignments: pd.DataFrame, max_misalignments: int
) -> set[int]:
    """Counts how often a primer misaligns and returns a set of primers that misalign more than max_misalignments times"""
    if max_misalignments < 0:
        return set()

    # Count the number of misalignments for each primer
    misalignments = alignments[
        alignments["primer_id"].map(alignments["primer_id"].value_counts())
        > max_misalignments
    ]

    # Return a set of primer_ids
    return set(misalignments["primer_id"].unique())


def __select_adjacent_primers(adjacent_primers_for_primer: dict[int, list]) -> set[int]:
    """Select misaligned primers which will be removed in the database"""

    # Transform adjacent_primers_for_primer into a list of tuples
    # Each tuple contains the primer_id and the list of adjacent primers
    connections = []
    for primer_id, adjacent_primers_infos in adjacent_primers_for_primer.items():
        primer_id = int(primer_id)
        for adjacent_primer_info in adjacent_primers_infos:
            adjacent_primer_id = int(adjacent_primer_info["primer_id"])
            connections.append((primer_id, adjacent_primer_id))

    # Filter out duplicate connections
    connections = list(set(connections))

    # Create an undirected graph from the connections
    graph = Graph(connections, directed=False)

    # Find the minimum number of nodes to remove to have no edges left -> vertices cover problem
    primer_ids_to_discard = graph.vertex_cover_approx()

    return set(primer_ids_to_discard)


def __mark_primers_as_discarded(db: DBHandler, primer_ids_to_discard: set[int]) -> None:
    """Marks the primers as discarded in the database"""

    # Create the query
    query = """
    UPDATE proto_primers
    SET discarded = 1
    WHERE id = ?
    """

    # Transform the set of primer_ids to discard into a list of tuples
    primer_ids_to_discard = [(primer_id,) for primer_id in primer_ids_to_discard]

    # Execute the query for each primer_id
    db.executemany(query, primer_ids_to_discard)


def __write_primers_as_discarded_to_file(
    primer_ids_to_discard: set[int], output: str
) -> None:
    """Writes the discarded primers to a file"""

    with open(output, "w") as f:
        for primer_id in primer_ids_to_discard:
            f.write(f"{primer_id}\n")


def main():
    logging.info("Starting evaluation of primers against target species")

    args = __get_args()
    db = DBHandler(args.db)

    alignments = __get_alignments(db, args)
    logging.info(f"Found {alignments.shape[0]} alignments to evaluate")

    adjacent_primers_for_primer = __calculate_adjacent_alignments(
        alignments, args.adjacency_limit
    )

    # Select primers to discard by unioning the sets of primers to discard
    primers_to_discard = set()
    primers_to_discard |= __select_adjacent_primers(adjacent_primers_for_primer)
    primers_to_discard |= __select_misaligned_primers(
        alignments, args.max_misalignments
    )

    logging.info(f"Marking {len(primers_to_discard)} primers as discarded")

    __mark_primers_as_discarded(db, primers_to_discard)
    __write_primers_as_discarded_to_file(primers_to_discard, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
