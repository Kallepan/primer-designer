import argparse
from db import DBHandler
import pandas as pd

MATCH_MULTIPLIER = 1
PRIMER_END_CUTOFF = 5
CLOSE_TO_3_PRIME_FACTOR = 0.5
FAR_FROM_3_PRIME_FACTOR = 0.1
BASE_PENALTY = 100
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score alignments from bowtie")
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
        help="Pool number to evaluate",
    )
    
    # optional arguments
    parser.add_argument(
        "--match_multiplier",
        type=int,
        default=MATCH_MULTIPLIER,
        help=f"Penalty multiplier for each alignment match not at the original position. Default: {MATCH_MULTIPLIER}",
    )
    parser.add_argument(
        "--primer_end_cutoff",
        type=int,
        default=PRIMER_END_CUTOFF,
        help=f"Number of bases from the 3' End of the primer where the close_to_3_prime_factor is applied. Default: {PRIMER_END_CUTOFF}",
    )
    parser.add_argument(
        "--close_to_3_prime_factor",
        type=float,
        default=CLOSE_TO_3_PRIME_FACTOR,
        help=f"Penalty for each mismatch within primer_end_cutoff bases from the 3' End of the primer. Default: {CLOSE_TO_3_PRIME_FACTOR}",
    )
    parser.add_argument(
        "--far_from_3_prime_factor",
        type=float,
        default=FAR_FROM_3_PRIME_FACTOR,
        help=f"Penalty for each mismatch outside primer_end_cutoff bases from the 3' End of the primer. Default: {FAR_FROM_3_PRIME_FACTOR}",
    )
    parser.add_argument(
        "--base_penalty",
        type=int,
        default=BASE_PENALTY,
        help=f"Default penalty applied for each misalignment of the primer. Default: {BASE_PENALTY}",
    )
    return parser.parse_args()

def get_alignments(db: DBHandler, args: argparse.Namespace) -> pd.DataFrame:
    # Extract all alignments along with the sequence of the primer
    alignments = db.select(
        """
        SELECT alignments.id, alignments.primer_id, alignments.aligned_to, alignments.matches, alignments.sequence, alignments.mismatches_descriptor, proto_primers.sequence AS primer_sequence, proto_primers.strand AS primer_strand
        FROM alignments
        LEFT JOIN proto_primers ON alignments.primer_id = proto_primers.id
        WHERE alignments.pool = ? AND (
            -- Filter out correct alignments so they don't get scored
            alignments.matches <> 1 OR
            alignments.mismatches_descriptor IS NOT NULL OR
            proto_primers.strand <> alignments.aligned_to
        )
        ORDER BY alignments.id ASC
        """, (args.pool,)
    )
    df = pd.DataFrame(
        alignments,
        columns=["id", "primer_id", "aligned_to", "matches", "sequence", "mismatches_descriptor", "primer_sequence", "primer_strand"]
    )
    return df

def score_alignments(alignments: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    def calc_score(row: pd.Series, match_multiplier: int, primer_end_cutoff: int, close_to_3_prime_factor: int, far_from_3_prime_factor: int, base_penalty: int) -> float:
        # Calculate the score for each alignment
        # The alignments calculated here should all be incorrect alignments
        mismatch_descriptor = row["mismatches_descriptor"]
        sequence = row["primer_sequence"]
        sequence_len = len(sequence)
        matches = row["matches"]
        # Filter correct alignments
        if mismatch_descriptor is None:
            return matches * match_multiplier + base_penalty
        
        # Try to find out where in the primer the mismatch is located
        mismatches = mismatch_descriptor.split(",")
        badness = matches * match_multiplier + base_penalty
        for mismatch in mismatches:
            mismatch_pos = int(mismatch.split(":")[0])
            distance_from_3_prime = sequence_len - mismatch_pos
            # The closer the mismatch is to the 3' end of the primer, the more unlikely this alignment is, so the lower the penalty
            badness -= (close_to_3_prime_factor if distance_from_3_prime < primer_end_cutoff else far_from_3_prime_factor) * distance_from_3_prime
        badness = max(badness, 0)
        return badness

    alignments["score"] = alignments.apply(calc_score, axis=1, args=(args.match_multiplier, args.primer_end_cutoff, args.close_to_3_prime_factor, args.far_from_3_prime_factor, args.base_penalty))

    return alignments

def to_db(alignments: pd.DataFrame, db: DBHandler, args: argparse.Namespace) -> None:
    # Write the scores to the database and to a file
    alignments.to_csv(args.output, sep="\t", index=False)
    db.executemany(
        """
        UPDATE alignments
        SET score = ?
        WHERE id = ?
        """, alignments[["score", "id"]].values.tolist()
    )

def main():
    args = get_args()
    db = DBHandler(args.db)
    alignments = get_alignments(db, args)
    scored_alignments = score_alignments(alignments, args)
    to_db(scored_alignments, db, args)

if __name__ == "__main__":
    print("Running score_alignments.py")
    main()