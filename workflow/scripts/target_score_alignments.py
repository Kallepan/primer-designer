import argparse
from db import DBHandler
import pandas as pd

MISMATCH_WEIGHT = 0.5
ALIGNMENT_WEIGHT = 0.5
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
    
    # Optional Args
    parser.add_argument(
        "--mismatch_weight",
        type=float,
        default=MISMATCH_WEIGHT,
        help=f"Weight applied to each mismatch of the alignment. Default: {MISMATCH_WEIGHT}",
    )
    parser.add_argument(
        "--alignment_weight",
        type=float,
        default=ALIGNMENT_WEIGHT,
        help=f"Weight applied to each alignment to the strand. Default: {ALIGNMENT_WEIGHT}",
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
    alignments, column_names = db.select(
        """
        SELECT alignments.id, alignments.primer_id, alignments.aligned_to, alignments.matches, alignments.sequence, alignments.mismatches_descriptor, proto_primers.sequence AS primer_sequence, proto_primers.strand AS primer_strand
        FROM alignments
        LEFT JOIN proto_primers ON alignments.primer_id = proto_primers.id
        WHERE alignments.pool = ? AND (
            -- Filter out correct alignments so they don't get scored, default: 0.0
            alignments.matches <> 1 OR
            alignments.mismatches_descriptor IS NOT NULL OR
            proto_primers.strand <> alignments.aligned_to
        )
        ORDER BY alignments.id ASC
        """, (args.pool,)
    )
    df = pd.DataFrame(
        alignments,
        columns=column_names,
    )
    return df

def score_alignments(alignments: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    def calc_score(row: pd.Series, base_penalty: int, alignment_weight: int, mismatch_weight: float) -> float:
        # Calculate the score for each alignment
        # The alignments calculated here should all be mismatched alignments
        mismatch_descriptor = row["mismatches_descriptor"]
        sequence = row["primer_sequence"]
        sequence_len = len(sequence)
        matches = row["matches"]
        # If no mismatches are present, alignments get the base penalty
        if mismatch_descriptor is None:
            return matches * base_penalty * alignment_weight

        # Try to find out where in the primer the mismatch is located
        mismatches = mismatch_descriptor.split(",")
        badness = matches * base_penalty * alignment_weight
        distance_sum = 1.0
        for mismatch in mismatches:
            mismatch_pos = int(mismatch.split(":")[0])
            # Calculate the distance from the 3' end of the primer by taking the length of the sequence and subtracting the mismatch position index + 1
            distance_from_3_prime = sequence_len - (mismatch_pos + 1)
            distance_sum += distance_from_3_prime
        # Calculate the final formula
        badness -= 1/distance_sum * len(mismatches) * mismatch_weight
        badness = max(badness, 0)
        return badness

    alignments["score"] = alignments.apply(calc_score, axis=1, args=(args.base_penalty, args.alignment_weight, args.mismatch_weight,))

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
    print("Running score_alignments.py")
    args = get_args()
    db = DBHandler(args.db)
    alignments = get_alignments(db, args)
    scored_alignments = score_alignments(alignments, args)
    to_db(scored_alignments, db, args)

if __name__ == "__main__":
    main()