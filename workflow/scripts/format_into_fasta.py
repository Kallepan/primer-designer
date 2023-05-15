import json
import argparse

import os
from io import TextIOWrapper

def __parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    return parser.parse_args()


def __load_json(path: str):
    with open(path, "r") as f:
        return json.load(f)


def __get_sequence_from_primer(primer: dict) -> str:
    return primer.get("primer_sequence", "")


def __extract_primer(
    primers: list[dict], 
    file: TextIOWrapper, 
    region_name: str, 
    amplicon_name: str, 
    strand: str
) -> None:
    if not primers:
        return None

    for i, primer in enumerate(primers):
        sequence = __get_sequence_from_primer(primer)
        name_str = f">{region_name}|{amplicon_name}|{strand}|{i}\n"
        file.write(name_str)
        file.write(sequence)
        file.write("\n")


def __write_fasta(data: dict, file: TextIOWrapper):
    for region in data["regions"]:
        region_name = region["region_name"]

        for amplicon in region["amplicons"]:
            amplicon_name = amplicon["amplicon_name"]

            # Extract the sequence frome the primer and write to file
            # Name Encoding: region_name|amplicon_name|(forward || reverse)|index
            __extract_primer(
                amplicon["forward_primers"], file, region_name, amplicon_name, "forward"
            )
            __extract_primer(
                amplicon["reverse_primers"], file, region_name, amplicon_name, "reverse"
            )


def main():
    args = __parse_args()
    data = __load_json(args.input)
    try:
        file = open(args.output, "w")
        __write_fasta(data, file)
    except Exception as e:
        print(f"Error writing to file {args.output}")
        print(e)
    finally:
        file.close()


if __name__ == "__main__":
    print("Formatting primers into fasta")
    main()
