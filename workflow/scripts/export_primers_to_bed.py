from Bio import SeqIO

import argparse
import logging
import json

logging.basicConfig(level=logging.INFO)


def __load_json(path: str) -> list[dict]:
    # Load json file
    with open(path, "r") as f:
        return json.load(f)


def __get_args() -> argparse.Namespace:
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Formats a set of primers into a bed file"
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input file",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Fasta file of the target genome",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )

    return parser.parse_args()


def __get_chromosome(fasta_file: str) -> str:
    """Gets the chromosome name from the fasta file"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        return record.id

def __format_into_bed(args: argparse.Namespace, chromosome: str, data: list[dict]) -> None:
    """Formats the input dictionary into a bed file"""
    with open(args.output, "w") as f:
        f.write("chrom\tstart\tend\tname\tscore\tstrand\n")

    # Simply loop over the data and write it to the file
    for pool in data:
        for primer_pair in pool["primer_pairs"]:
            forward_primer = primer_pair["forward_primer"]
            reverse_primer = primer_pair["reverse_primer"]
            name = primer_pair["amplicon_name"]

            # if the primer could not be aligned, skip it
            # We set the position to -1 if it could not be aligned
            if forward_primer["position"] < 0:
                logging.warning(
                    f"Forward primer: {forward_primer['id']} could not be aligned."
                )
                continue
            if reverse_primer["position"] < 0:
                logging.warning(
                    f"Reverse primer: {reverse_primer['id']} could not be aligned."
                )
                continue

            start = forward_primer["position"]
            end = reverse_primer["position"] + reverse_primer["length"]
            with open(args.output, "a") as f:
                f.write(f"{chromosome}\t{start}\t{end}\t{name}\t0\t+\n")


def main():
    logging.info("Starting script")
    args = __get_args()
    data = __load_json(args.input)
    
    # Get the chromosome name from the fasta file
    chromosome = __get_chromosome(args.fasta)

    __format_into_bed(args, chromosome, data)

    logging.info("Script finished")
    logging.info(
        "To view the primers in IGV, load in the original fasta file of the target genome and the bed file generated by this script."
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(e)
        raise e
