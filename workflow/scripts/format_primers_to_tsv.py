"""
Formats the primers into an output tsv file:
All sequences are outputted in 5' -> 3'
    - amplicon_name: name of the amplicon
    - forward_sequence: sequence of the forward primer
    - reverse_sequence: sequence of the reverse primer
"""

import logging
import argparse
import json

import pandas as pd

from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Formatting primers into tsv")

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to final primers json file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output file",
    )

    return parser.parse_args()


def __format_to_tsv(input_file: str, output_file: str) -> None:
    with open(input_file, "r") as f:
        primers = json.load(f)

    # Iterate through each primer set and merge into one dataframe
    logging.info("All sequences are outputted in 5' -> 3'")
    column_names = [
        "amplicon_name",
        "forward_primer_sequence",
        "reverse_primer_sequence",
    ]
    data = []
    for primer_set in primers:
        primer_pairs = primer_set.get("primer_pairs", [])

        # Iterate through each primer pair
        for primer_pair in primer_pairs:
            amplicon_name = primer_pair.get("amplicon_name", "")
            forward_primer_sequence = primer_pair.get("forward_primer", {}).get(
                "sequence", ""
            )
            reverse_primer_sequence = primer_pair.get("reverse_primer", {}).get(
                "sequence", ""
            )

            data.append(
                [amplicon_name, forward_primer_sequence, reverse_primer_sequence]
            )

    df = pd.DataFrame(data, columns=column_names)

    # write to tsv
    logging.info(f"Writing {df.shape[0]} rows to {output_file}")
    df.to_csv(output_file, index=False, sep="\t")


def main():
    logging.info("Formatting primers into tsv")

    args = __get_parser()

    __format_to_tsv(args.input, args.output)

    logging.info("Done")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        raise e
