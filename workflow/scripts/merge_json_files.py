import argparse
import json
import logging

logging.basicConfig(level=logging.INFO)


def __get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merges json files with object as root"
    )
    parser.add_argument(
        "--input", type=str, help="Input files", nargs="+", required=True
    )
    parser.add_argument("--output", type=str, help="Output file", required=True)
    return parser.parse_args()


def __merge_json_files(input_files: list, output_file: str) -> None:
    # Merge json output files
    output = []

    # Load all input files
    for input_file in input_files:
        with open(input_file, "r") as f:
            output.append(json.load(f))

    # Write output file
    with open(output_file, "w") as f:
        json.dump(output, f, indent=4)


def main():
    logging.info("Merging saddle output files")

    args = __get_parser()
    __merge_json_files(args.input, args.output)

    logging.info("Done")


if __name__ == "__main__":
    main()
