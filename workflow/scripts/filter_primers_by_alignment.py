import argparse
import pandas as pd
import json
import sys

DEFAULT_ADJACENCY_LIMIT = 500

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter primers by alignment")
    parser.add_argument(
        "--alignment",
        type=str,
        required=True,
        help="Alignment file with primer sequences",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file with primer filtered sequences",
    )
    parser.add_argument(
        "--primers",
        type=str,
        required=True,
        help="JSON file with original primer sequences",
    )
    parser.add_argument(
        "--adjacency_limit",
        type=int,
        default=DEFAULT_ADJACENCY_LIMIT,
        help=f"Limit within which primers of opposite strand are considered adjancent enough to be problematic. Default: {DEFAULT_ADJACENCY_LIMIT}",
    )

    return parser.parse_args()


def __write_output(primers: dict, output: str) -> None:
    with open(output, "w") as file:
        json.dump(primers, file, indent=4)


def __get_primers_to_remove(alignment: pd.DataFrame, args: argparse.Namespace) -> list:
    """
    Takes in a dataframe containing the parsed output from bowtie
    and returns a list of primers to remove

    Algorithm:
        1. Find all primers not aligning to the original site:
            - Primers aligning multiple times
            - Primers with mismatches
            - Primers aligning to the wrong strand
        2. For each primer, find all primers aligning to the other strand within a certain distance
        3. Calculate a badness score for each primer considering:
            - the amount of mismatches
            - misalignments for the primer
            - the amount of adjacent primers (inverse the distance)
        4. Append the score in the filtered_proto_primers.json file and output
    """
    problematic_primers = alignment.loc[
        (alignment["matches"] > 1) |
        (alignment["strand"] != alignment["primer_strand"]) |
        (alignment["mismatches_descriptor"].isna() == False)
    ]
    problematic_primers.reset_index(inplace=True)

     # check each misaligned_primer for adjacent primers -> opposite strand, within adjacency_limit
    print(f"Found {len(problematic_primers)} misaligned primers")
    for _, problematic_primer in problematic_primers.iterrows():
        # TODO: Implement statistics logging for end-user
        adjacency_limit = args.adjacency_limit if problematic_primer["primer_strand"] == "forward" else -1*args.adjacency_limit
        opposite_strand = "reverse" if problematic_primer["primer_strand"] == "forward" else "forward"
        adjacent_primers = alignment.loc[
            (alignment["primer_amplicon"] != problematic_primer["primer_amplicon"]) &
            (alignment["strand"] == opposite_strand) &
            (alignment["position"] <= problematic_primer["position"] + adjacency_limit) &
            (alignment["position"] >= problematic_primer["position"])
        ]
        print(f"Found {len(adjacent_primers)} adjacent primers for primer {problematic_primer['primer']}")
        if len(adjacent_primers) < 0:
            continue


def main():
    args = get_args()

    # Read alignment file
    # Columns: primer	strand	chromosome	position	sequence	read_quality	matches	primer_region	primer_amplicon	primer_strand	primer_id
    dtypes = {
        "primer": "str",
        "strand": "str",
        "chromosome": "str",
        "position": "int",
        "sequence": "str",
        "read_quality": "str",
        "matches": "int",
        "primer_region": "str",
        "primer_amplicon": "str",
        "primer_strand": "str",
        "primer_id": "str",
    }
    alignment = pd.read_csv(args.alignment, sep="\t", header=0, dtype=dtypes)

    # Read primer file
    with open(args.primers, "r") as file:
        original_primers = json.load(file)

    # Filter primers by alignment

    # Write output
    __get_primers_to_remove(alignment, args)
    __write_output(original_primers, args.output)


if __name__ == "__main__":
    print("Starting primer filtering")
    main()
    try:
        pass
    except Exception as e:
        print(e)
        sys.exit(1)
