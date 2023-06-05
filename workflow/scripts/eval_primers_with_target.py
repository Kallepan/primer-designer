import argparse
import sys
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
        "--hard_filter",
        action="store_true",
        help=f"Enabling this option removes primers which misalign and have primers on the reverse complement strand adjacent to the alignment in the direction of amplification."
    )

    return parser.parse_args()

def calculate_badness(db: DBHandler, args: argparse.Namespace) -> None:
    """
    Calculate the self-badness for each proto_primer in the database and output the result to the database as well as a json file.
    The self badness is calculated using following algorithm:

        1. Find all primers not aligning to the original site:
            - Primers aligning multiple times
            - Primers with mismatches
            - Primers aligning to the wrong strand
        2. For each primer, find all primers aligning to the other strand within a certain distance not of the same amplicon
        3. Calculate a badness score for each primer considering:
            - the amount of mismatches
            - misalignments for the primer
            - the amount of adjacent primers (inverse the distance)
        4. Append the score to the database and json file
    """
    # Find all primers not aligning to the original site and associate them with primers aligning to the other strand within a certain distance
    # 
    problematic_primers = db.select(
        """
        -- Select all problematic primers
        WITH problematic_primers AS (
            SELECT alignments.id AS problematic_id, alignments.pool AS problematic_pool, proto_primers.amplicon_name AS problematic_amplicon_name, alignments.position AS problematic_position, alignments.matches AS problematic_matches, alignments.mismatches_descriptor AS problematic_mismatches_descriptor, alignments.aligned_to AS problematic_aligned_to, proto_primers.strand AS problematic_original_strand
            FROM alignments 
            LEFT JOIN proto_primers ON alignments.id = proto_primers.id
            WHERE alignments.pool = ? AND (
                alignments.matches > 1 OR 
                alignments.mismatches_descriptor IS NOT NULL OR 
                alignments.aligned_to <> proto_primers.strand
            )
            ORDER BY alignments.id ASC
        ), -- now add the adjacent primers to these primers
        formatted_alignments AS (
            SELECT alignments.id, alignments.pool, proto_primers.amplicon_name, alignments.position, alignments.matches, alignments.mismatches_descriptor, alignments.aligned_to AS alignment_aligned_to, proto_primers.strand AS original_strand
            FROM alignments
            LEFT JOIN proto_primers ON alignments.id = proto_primers.id
            WHERE alignments.pool = ?
        ) -- now we have the alignments with their corresponding strand & amplicon_name
        SELECT * 
        FROM problematic_primers
        LEFT JOIN formatted_alignments AS adjacent_primers ON
        problematic_amplicon_name <> adjacent_primers.amplicon_name AND
        problematic_aligned_to <> adjacent_primers.alignment_aligned_to AND
        problematic_pool = adjacent_primers.pool AND
        -- The following code should be done in python not in sql, because it is not very readable
        CASE
            WHEN problematic_aligned_to = 'forward' THEN
                adjacent_primers.position >= problematic_position AND
                adjacent_primers.position <= problematic_position + ? -- adjacency_limit
            WHEN problematic_aligned_to = 'reverse' THEN
                adjacent_primers.position <= problematic_position AND
                adjacent_primers.position >= problematic_position - ? -- adjacency_limit
        END
        """, (args.pool, args.pool, args.adjacency_limit, args.adjacency_limit)
    )

    problematic_primers_df = pd.DataFrame(
        problematic_primers,
        columns = [
            "problematic_id",
            "problematic_pool",
            "problematic_amplicon_name",
            "problematic_position",
            "problematic_matches",
            "problematic_mismatches_descriptor",
            "problematic_alignment_strand",
            "problematic_original_strand",
            "adjacent_id",
            "adjacent_pool",
            "adjacent_amplicon_name",
            "adjacent_position",
            "adjacent_matches",
            "adjacent_mismatches_descriptor",
            "adjacent_alignment_strand",
            "adjacent_original_strand",
        ]
    )
    # TODO: DEBUG
    problematic_primers_df.to_csv("tmp/problematic_primers.tsv", sep="\t", index=False)

    # Now calculate the score of each problematic primer and append it to the database
    # Take into consideration if any adjacent primers are present

    grouped_primers = problematic_primers_df.groupby("problematic_id")
    summary_df = grouped_primers.agg(
        id=("problematic_id", "first"),
        adjacent_alignments=("adjacent_id", "count"), # count the amount of adjacent primers
        total_problematic=("problematic_id", "count"), # count the amount of alignments + adjacent primers 
    )
    # Export summary scores to dataframe and csv
    summary_df.to_csv(args.output , sep="\t", index=False)

    # Export the summary df to the database    
    # Update the proto_primers table with the badness score
    # Drop the temporary table if it exists

def main():
    args = get_args()
    db = DBHandler(args.db)

    calculate_badness(db, args)

if __name__ == "__main__":
    print("Evaluating primers against target species")
    main()
