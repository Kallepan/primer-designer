import logging
import argparse
import sys
import pandas as pd

from collections import defaultdict
from handlers import DBHandler
from handlers import Graph

logging.basicConfig(level=logging.INFO)

DEFAULT_ADJACENCY_LIMIT = 500
DEFAULT_BASES_TO_IGNORE = 3
DEFAULT_MAX_MISMATCHES = 2


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
        help=f"Number of bases at the 3' end of the primer to ignore when evaluating alignments. Default: {DEFAULT_BASES_TO_IGNORE}",
    )
    parser.add_argument(
        "--max_mismatches",
        type=int,
        default=DEFAULT_MAX_MISMATCHES,
        help=f"Maximum number of mismatches allowed for an alignment to be considered during calculation. Default: {DEFAULT_MAX_MISMATCHES}",
    )

    return parser.parse_args()


def __get_alignments_from_species(
    args: argparse.Namespace, db: DBHandler, pool: str, species: str
) -> pd.DataFrame:
    """
    Fetch all alignments from a given species and pool. At first I wanted to do an sql join but it was too slow. I suppose joining 1_000_000 with 1_000_000 rows is too much for sqlite.
    I could maybe groupby the chromosome and reduce the amount of rows to join to 50_000. But I think it's better to just fetch all the alignments and then filter them in python. See the __filter_alignments function.
    """

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
        alignments.pool = ? AND 
        alignments.species = ? AND 
        ((
        -- do not filter out alignments with mismatches_descriptor if the primer is too short or if mismatches_descriptor is NULL
            mismatches_descriptor IS NULL OR 
            LENGTH(proto_primers.sequence) < 10
        ) OR (
        -- filter out alignments with mismatches_descriptor if the primer is long enough and mismatches_descriptor is not NULL
    """

    # Add the bases to ignore to the query
    for i in range(1, args.bases_to_ignore + 1):
        query += f""" INSTR(alignments.mismatches_descriptor, CAST(LENGTH(proto_primers.sequence) - {i} AS TEXT)) = 0 AND """

    # Finish the query
    query += """
    1 = 1))
    ORDER BY alignments.chromosome ASC, alignments.position ASC
    """
    alignments, columns = db.select(
        query,
        (pool, species),
    )

    return pd.DataFrame(alignments, columns=columns)


def __calculate_adjacent_alignments(
    alignments: pd.DataFrame, adjacency_limit: int, max_mismatches: int
) -> dict[str, list[int]]:
    """
    For each alignment, find the next alignment that is within the adjacency limit.
    Keep track of the adjacent alignments for each primer in a dict.
    """

    # Split alignments into tables for each chromosome
    alignments_by_chromosome = alignments.groupby("chromosome")

    adjacent_primers_for_primer: dict[str, list[int]] = defaultdict(list)
    for chromosome, alignments in alignments_by_chromosome:
        logging.info(
            f"Processing chromosome {chromosome} with {alignments.shape[0]} alignments"
        )

        # Sort alignments by position position
        alignments = alignments.sort_values("position")

        # Filter the alignments to only keep the ones that have no mismatches descriptor or have a mismatches descriptor with less than max_mismatches mismatches
        # we can simply count the occurences of ">" in the mismatches descriptor to know how many mismatches there are
        alignments = alignments[
            (alignments["mismatches_descriptor"].isnull())
            | (alignments["mismatches_descriptor"].str.count(">") <= max_mismatches)
        ]

        # Split alignments into forward and reverse alignments
        forward_alignments = alignments[alignments["primer_strand"] == "forward"]
        reverse_alignments = alignments[alignments["primer_strand"] == "reverse"]

        # For each alignment, in the dataframes, find all alignments of the complementary strand that are within the adjacency limit
        for _, alignment in forward_alignments.iterrows():
            # Filter out complementary alignments that are not within the adjacency limit
            adjacent_alignments = reverse_alignments[
                (reverse_alignments["position"] >= alignment["position"])
                & (
                    reverse_alignments["position"]
                    <= alignment["position"] + adjacency_limit
                )
            ]

            # if adjacent_alignments is empty, we can continue
            if adjacent_alignments.shape[0] == 0:
                continue

            # Extract the information we need from the adjacent_alignments
            adjacent_alignments_info = adjacent_alignments["primer_id"].tolist()

            # Add the adjacent alignments to the dict
            alignment_id = int(alignment["primer_id"])
            adjacent_primers_for_primer[alignment_id].extend(adjacent_alignments_info)

        # Repeat the same process for the reverse alignments
        for _, alignment in reverse_alignments.iterrows():
            # Filter out complementary alignments that are not within the adjacency limit
            adjacent_alignments = forward_alignments[
                (
                    forward_alignments["position"]
                    >= alignment["position"] - adjacency_limit
                )
                & (forward_alignments["position"] <= alignment["position"])
            ]

            # if adjacent_alignments is empty, we can continue
            if adjacent_alignments.shape[0] == 0:
                continue

            # Extract the information we need from the adjacent_alignments
            adjacent_alignments_info = adjacent_alignments["id"].tolist()

            # Add the adjacent alignments to the dict
            alignment_id = int(alignment["id"])
            adjacent_primers_for_primer[alignment_id].extend(adjacent_alignments_info)

    return adjacent_primers_for_primer


def __select_primers_to_exlude(
    adjacent_primers_for_primer: dict[int, list[int]]
) -> list[int]:
    """Determines which primers to mark as dicarded in the databases based on their adjacency"""

    # Transform adjacent_primers_for_primer into a list of tuples
    # Each tuple contains the primer_id and the list of adjacent primers
    connections = []
    for primer_id, adjacent_primer_ids in adjacent_primers_for_primer.items():
        primer_id = int(primer_id)
        for adjacent_primer_id in adjacent_primer_ids:
            connections.append((primer_id, adjacent_primer_id))

    # Filter out duplicates
    connections = list(set(connections))

    # Create an undirected graph from the connections
    graph = Graph(connections, directed=False)

    # Find the minimum number of nodes to remove to have no edges left -> vertices cover problem
    primer_ids_to_discard = graph.vertex_cover_approx()

    return primer_ids_to_discard


def __mark_primers_as_discarded(
    db: DBHandler, primer_ids_to_discard: list[int]
) -> None:
    """Marks the primers as discarded in the database"""

    # Create the query
    query = """
    UPDATE proto_primers
    SET discarded = 1
    WHERE id = ?
    """

    primer_ids_to_discard = [(primer_id,) for primer_id in primer_ids_to_discard]

    # Execute the query for each primer_id
    db.executemany(query, primer_ids_to_discard)


def main():
    logging.info("Scoring primers against foreign genomes")
    args = __get_args()

    db = DBHandler(args.db)

    alignments = __get_alignments_from_species(args, db, args.pool, args.species)

    logging.info(f"Found {alignments.shape[0]} alignments")
    adjacent_primers_for_primer = __calculate_adjacent_alignments(
        alignments, args.adjacency_limit, args.max_mismatches
    )
    primers_to_discard = __select_primers_to_exlude(adjacent_primers_for_primer)

    logging.info(f"Marking {len(primers_to_discard)} primers as discarded")
    __mark_primers_as_discarded(db, primers_to_discard)

    import json

    with open(args.output, "w") as f:
        json.dump(adjacent_primers_for_primer, f, indent=4)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
