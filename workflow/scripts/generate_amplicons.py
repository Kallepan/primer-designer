"""
Module to load in fasta file of an organism and extract the regions of interest (loci) from the genome.
To use: python load_fasta.py <path_to_fasta_file> <path_to_loci_file>
"""

import sys
import os
import argparse

from collections import defaultdict

from Bio import SeqIO
import pandas as pd


def _extract_seq_record_from_fasta_file(path_to_fasta_file: str) -> SeqIO.SeqRecord:
    """Read a fasta file and return a list of sequences"""

    with open(path_to_fasta_file, "r") as handle:
        fa = list(SeqIO.parse(handle, "fasta"))

    if len(fa) == 0:
        raise Exception("No sequences found in fasta file")

    if len(fa) > 1:
        raise Exception("More than one sequence found in fasta file")

    return fa[0]


def _load_loci(path_to_loci_file: str) -> pd.DataFrame:
    """Read a loci csv file and return a list of loci"""

    df = pd.read_csv(
        path_to_loci_file,
        sep=",",
        header=0,
        names=["loci", "start", "end"],
        dtype={"loci": str, "start": int, "end": int},
    )

    return df


def _extract_region(seq, start: int, end: int, buffer: int = 100) -> SeqIO.SeqRecord:
    """Extract a region of genome from a sequence"""

    # TODO: What to do when the region is smaller than the amplicon?

    if start > end:
        return seq[end - buffer : start + buffer]

    return seq[start - buffer : end + buffer]


def _split_region_into_amplicons(
    region, amplicon_size: int = 1000, amplicon_buffer: int = 100
) -> list[SeqIO.SeqRecord]:
    """Split a region of the genome into amplicons"""
    amplicon_size = amplicon_size - amplicon_buffer
    amplicons = []

    # process initial amplicon separately
    amplicons.append(region[0 : amplicon_size + amplicon_buffer])

    for i in range(amplicon_size, len(region) - amplicon_size, amplicon_size):
        amplicons.append(region[i - amplicon_buffer : i + amplicon_size])

    # process final amplicon separately
    amplicons.append(region[-amplicon_size - amplicon_buffer :])

    return amplicons


def _write_amplicons_to_file(amplicons: dict, path_to_output_file: str) -> None:
    """Write amplicons to file"""
    with open(path_to_output_file, "w") as handle:
        for region_name, region_amplicons in amplicons.items():
            for i, amplicon in enumerate(region_amplicons):
                handle.write(f">{region_name}_{i}")
                handle.write("\n")
                handle.write(f"{amplicon}")
                handle.write("\n")


def get_parser() -> argparse.ArgumentParser:
    """
    Get parser for script. Is needed if the script is called from the command line.
    """

    parser = argparse.ArgumentParser(
        prog="A script to load in fasta file of an organism and extract the regions of interest (loci) from the genome and split these regions into amplicons of desired length.",
        description="Load in fasta file of an organism and extract the regions of interest (loci) from the genome.",
    )

    parser.add_argument("fasta_file", type=str, help="Path to fasta file")
    parser.add_argument(
        "loci_file",
        type=str,
        help="Path to loci file. Format: CSV. Columns: loci, start, end",
    )
    parser.add_argument(
        "output_file", type=str, help="Path to output file in fasta format"
    )
    parser.add_argument(
        "-s",
        "--amplicon_size",
        type=int,
        help="Size of amplicons. Default: 1000.",
        required=False,
    )
    parser.add_argument(
        "-b",
        "--amplicon_buffer",
        type=int,
        help="Buffer size for amplicons. Default: 100",
        required=False,
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    fasta_file = args.fasta_file
    loci_file = args.loci_file
    output_file_path = args.output_file
    amplicon_size = args.amplicon_size
    amplicon_buffer = args.amplicon_buffer

    # Check if fasta file exists
    if not os.path.exists(fasta_file):
        raise Exception("Fasta file does not exist")

    # Check if loci file exists
    if not os.path.exists(loci_file):
        raise Exception("Loci file does not exist")

    if os.path.exists(output_file_path):
        raise Exception("Output file already exists")

    if not amplicon_buffer:
        amplicon_buffer = 100

    if not amplicon_size:
        amplicon_size = 1000

    if amplicon_size < amplicon_buffer:
        raise Exception("Amplicon size must be greater than amplicon buffer")

    # Load fasta file
    seq_record = _extract_seq_record_from_fasta_file(fasta_file)

    # Load loci file
    loci = _load_loci(loci_file)

    # Extract regions of interest
    regions = defaultdict(None)
    for i, row in loci.iterrows():
        loci_name = row["loci"]
        start = row["start"]
        end = row["end"]

        regions[loci_name] = _extract_region(seq_record.seq, start, end)

        # print message if region is smaller than amplicon
        if len(regions[loci_name]) < amplicon_size:
            print(
                f"Warning: Region {loci_name} is smaller than amplicon size of by {amplicon_size - len(regions[loci_name])} bp."
            )

    # Split regions into amplicons
    amplicons = defaultdict(list)
    for region_name, region in regions.items():
        amplicons[region_name] = _split_region_into_amplicons(
            region, amplicon_size, amplicon_buffer
        )

    # Write amplicons to file
    _write_amplicons_to_file(amplicons, output_file_path)

    return seq_record, amplicons


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"Error: {e}")
        sys.exit(1)
