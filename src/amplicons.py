"""
Module to load in fasta file of an organism and extract the regions of interest (loci) from the genome.
To use: python load_fasta.py <path_to_fasta_file> <path_to_loci_file>
"""

import sys
import os
import utils
import argparse

from collections import defaultdict

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="A script to load in fasta file of an organism and extract the regions of interest (loci) from the genome and split these regions into amplicons of desired length.",
        description="Load in fasta file of an organism and extract the regions of interest (loci) from the genome.",
    )

    parser.add_argument('fasta_file', type=str, help='Path to fasta file')
    parser.add_argument('loci_file', type=str, help='Path to loci file. Format: CSV. Columns: loci, start, end')
    parser.add_argument('output_file', type=str, help='Path to output file in fasta format')
    parser.add_argument('-s', '--amplicon_size', type=int, help='Size of amplicons. Default: 1000', required=False)
    parser.add_argument('-b', '--amplicon_buffer', type=int, help='Buffer size for amplicons. Default: 50', required=False)
    
    return parser

def format_fasta_into_amplicons():
    parser = get_parser()
    args = parser.parse_args()

    fasta_file = args.fasta_file
    loci_file = args.loci_file
    output_file_path = args.output_file
    amplicon_size = args.amplicon_size
    amplicon_buffer = args.amplicon_buffer

    # Check if fasta file exists
    if not os.path.exists(fasta_file):
        print("Fasta file does not exist")
        sys.exit(1)

    # Check if loci file exists
    if not os.path.exists(loci_file):
        print("Loci file does not exist")
        sys.exit(1)

    if os.path.exists(output_file_path):
        print("Output file already exists")
        sys.exit(1)

    if not amplicon_buffer:
        amplicon_buffer = 50

    # Load fasta file
    seq_record = utils.extract_seq_record_from_fasta_file(fasta_file)

    # Load loci file
    loci = utils.load_loci(loci_file)

    # Extract regions of interest
    regions = defaultdict(None)
    for i, row in loci.iterrows():
        loci_name = row["loci"]
        start = row["start"]
        end = row["end"]

        regions[loci_name] = utils.extract_region(seq_record.seq, start, end)

    # Split regions into amplicons
    amplicons = defaultdict(list)
    for region_name, region in regions.items():
        amplicons[region_name] = utils.split_region_into_amplicons(region, amplicon_size, amplicon_buffer)

    # Write amplicons to file
    utils.write_amplicons_to_file(amplicons, output_file_path)

    return seq_record, amplicons

if __name__ == "__main__": # pragma: no cover
    format_fasta_into_amplicons()