"""
Note that this script is not used in the workflow, but is provided as an example
"""

import argparse
import os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Split amplicons into separate files")
parser.add_argument("input", help="Input file")
parser.add_argument("output", help="Output directory")

args = parser.parse_args()

if not os.path.exists(args.input):
    sys.stderr.write(f"Input file {args.input} does not exist")
    exit(1)

if not os.path.exists(args.output):
    os.mkdir(args.output)

with open(args.input, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        with open(os.path.join(args.output, record.id + ".fasta"), "w") as out:
            SeqIO.write(record, out, "fasta")
