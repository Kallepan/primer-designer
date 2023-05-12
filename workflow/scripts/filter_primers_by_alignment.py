import argparse
import pandas as pd
import json
import sys


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
        help="File with original primer sequences",
    )

    return parser.parse_args()


def __write_output(primers: dict, output: str) -> None:
    with open(output, "w") as file:
        json.dump(primers, file, indent=4)


def __get_primers_to_remove(alignment: pd.DataFrame) -> list:
    # TODO implement method to get primers to remove
    return []


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
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
