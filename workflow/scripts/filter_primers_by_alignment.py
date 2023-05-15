import argparse
import pandas as pd
import json
import sys

DEFAULT_MIN_PRIMER_DISTANCE = 500

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
        "--min_primer_distance",
        type=int,
        default=DEFAULT_MIN_PRIMER_DISTANCE,
        help=f"Distance in bp to consider two primers as adjacent. Default: {DEFAULT_MIN_PRIMER_DISTANCE}",
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
        1. For each primer, check if it aligns more than once
        2. If yes, check if the opposite strand contains other close primers
        3. If no, continue, else check if it is the only primer of this amplicon in the dataframe
        4. If yes, add the reverse primers to the list of primers to remove, else add to list of primers to remove
    """
    min_primer_distance = args.min_primer_distance
    multiple_alignments = alignment[alignment["matches"] > 1]

    primers_to_remove = []
    for _, row in multiple_alignments.iterrows():
        pass

    return primers_to_remove


def main():
    args = get_args()

    # Read alignment file
    alignment = pd.read_csv(args.alignment, sep="\t", header=0)

    # Read primer file
    with open(args.primers, "r") as file:
        original_primers = json.load(file)

    # Filter primers by alignment

    # Write output
    # TODO change method to write output
    __write_output(original_primers, args.output)


if __name__ == "__main__":
    print("Starting primer filtering")
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
