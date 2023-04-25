import os
import sys
import json
import argparse
import yaml

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd

from handler import PrimerGenerator

DEFAULT_AMPLICON_SIZE=200
DEFAULT_MIN_AMPLICON_OVERLAP=50
DEFAULT_MAX_AMPLICON_OVERLAP=150
DEFAULT_OVERLAP_INCREASE=10
DEFAULT_OVERLAP_SHIFT=10

def __load_primer3_settings(path: str) -> str:
    """
    Load the primer3 settings from the config file
    """
    KEYS_TO_IGNORE = ["PRIMER_PRODUCT_SIZE_RANGE", "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"]
    settings = {}

    with open(path, "r") as handle:
        settings = yaml.safe_load(handle)
    
    settings_str = ""

    for k in KEYS_TO_IGNORE:
        settings.pop(k, None)

    for k, v in settings.items():
        settings_str += f"{k}={v}\n"

    settings_str += "="

    return settings_str 

def __load_regions(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=",",
        header=0,
        names=["loci", "start", "end"],
        dtype={"loci": str, "start": int, "end": int},
    )

    return df

def __extract_sequence_record_from_fasta(path: str) -> SeqRecord:
    with open(path, "r") as file:
        seq = list(SeqIO.parse(file, "fasta"))
    
    if len(seq) == 0:
        raise Exception("No sequences found in fasta file")
    
    if len(seq) > 1:
        raise Exception("More than one sequence found in fasta file")
    
    return seq[0]

def get_parser() -> argparse.ArgumentParser:
    """
        Returns the parser for the script.
    """

    parser = argparse.ArgumentParser(
        description=
        """
        Generate a list of amplicons and generate Proto Primers for each amplicon from a sequence in a fasta file.
        Save the list of amplicons to a fasta file as well as the Proto Primers to a json.
        """
    )
    parser.add_argument(
        "-f", 
        "--fasta",
        type=str,
        required=True,
        help="Path to the fasta file containing the sequence to be used to generate the amplicons."
    )
    parser.add_argument(
        "-c",
        "--config_file",
        type=str,
        required=True,
        help="Path to the primer3 config file."
    )
    parser.add_argument(
        "-r",
        "--regions",
        type=str,
        required=True,
        help="Path to the csv file containing the regions to be used to generate the amplicons."
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Path to the output file to be used to save the amplicons. Use the .json format."
    )
    parser.add_argument(
        "-t",
        "--temp_dir",
        type=str,
        help="Path to the temporary directory to be used when generating the amplicons."
    )
    parser.add_argument(
        "--amplicon_size",
        type=int,
        default=DEFAULT_AMPLICON_SIZE,
        help=f"The size of the amplicons to be generated. Default is {DEFAULT_AMPLICON_SIZE}."
    )
    parser.add_argument(
        "--min_overlap",
        type=int,
        default=DEFAULT_MIN_AMPLICON_OVERLAP,
        help=f"The minimum overlap to be used when generating the amplicons. Default is {DEFAULT_MIN_AMPLICON_OVERLAP}."
    )
    parser.add_argument(
        "--max_overlap",
        type=int,
        default=DEFAULT_MAX_AMPLICON_OVERLAP,
        help=f"The maximum overlap to be used when generating the amplicons. Default is {DEFAULT_MAX_AMPLICON_OVERLAP}."
    )
    parser.add_argument(
        "--overlap_increase",
        type=int,
        default=DEFAULT_OVERLAP_INCREASE,
        help="The amount by which the overlap will be increased when generating the amplicons and if the previous overlap returned no primer pairs."
    )
    parser.add_argument(
        "--overlap_shift",
        type=int,
        default=DEFAULT_OVERLAP_SHIFT,
        help="The amount by which the overlap will be shifted when generating the amplicons and if the previous overlap returned no primer pairs."
    )
    return parser

def __remove_duplicate_primers(list_of_primers: list[dict]) -> list:
    seen_primer = set()
    primers = []
    for primer_data in list_of_primers:
        if primer_data["sequence"] in seen_primer:
            continue
        seen_primer.add(primer_data["sequence"])
        primers.append(primer_data)
    return primers

def __generate_amplicon_and_primers_by_regions(
        region_name: str,
        region_start: int,
        region_end: int,
        sequence: Seq,
        amplicon_size: int,
        amplicon_max_overlap: int,
        amplicon_min_overlap: int,
        overlap_increase_increment: int,
        overlap_shift_increment: int,
        primer3_settings_path: str,
        temp_dir: str
    ) -> list:
    """
        Split a given region into amplicons by including a certain overlap in both directions.
    """
    primer3_settings = __load_primer3_settings(primer3_settings_path)

    amplicons = []

    amplicon_size = amplicon_size - (2 * amplicon_min_overlap)
    amplicon_offset = - amplicon_min_overlap
    index = 0
    while amplicon_offset < (region_end - region_start + amplicon_min_overlap):
        # Executed for each amplicon
        amplicon_overlap = amplicon_min_overlap
        amplicon_shift = 0
        primer_generation_iteration_index = 0
        amplicon_forward_primers = []
        amplicon_reverse_primers = []

        while True:
            # Shifting strategy primer generation
            primer_generation_iteration_index += 1
            amplicon_start  = amplicon_offset + amplicon_shift - amplicon_overlap
            amplicon_end    = amplicon_offset + amplicon_shift + amplicon_overlap + amplicon_size  

            if amplicon_shift >= amplicon_overlap:
                break
            
            amplicon_sequence = sequence[amplicon_start + region_start: amplicon_end + region_start]

            primer_generator = PrimerGenerator(
                region_name=region_name,
                amplicon_index=index,
                amplicon_sequence=amplicon_sequence,
                primer3_settings=primer3_settings,
                temp_dir=temp_dir,
                amplicon_overlap=amplicon_overlap,
                iteration_index=primer_generation_iteration_index
            )

            forward_primers, reverse_primers = primer_generator.run()
            # Append non-empty forward_primers and reverse_primers to the amplicons dict
            # if either is not empty and the list does not have entries break this loop
            
            amplicon_forward_primers.extend(forward_primers)
            amplicon_reverse_primers.extend(reverse_primers)

            if forward_primers and reverse_primers:
                break

            amplicon_shift += overlap_shift_increment
            continue

        # TODO: Add additional stragies for when no primers are found

        if not amplicon_forward_primers and not amplicon_reverse_primers:
            print(f"Failed to find primers for {region_name}-{index} in region {amplicon_offset}-{amplicon_offset + amplicon_size + amplicon_overlap}")

        amplicon_forward_primers = __remove_duplicate_primers(amplicon_forward_primers)
        amplicon_reverse_primers = __remove_duplicate_primers(amplicon_reverse_primers)

        amplicon = {
            "name": f"{region_name}-{index}",
            "forward_primers": amplicon_forward_primers,
            "reverse_primers": amplicon_reverse_primers,
        }
        amplicons.append(amplicon)
        amplicon_offset += amplicon_size
        index += 1

    return amplicons

def main():
    args = get_parser().parse_args()

    # Error Handling
    if not os.path.exists(args.fasta):
        raise Exception(f"The fasta file {args.fasta} does not exist.")
    if not os.path.exists(args.regions):
        raise Exception(f"The regions file {args.regions} does not exist.")
    if not os.path.exists(args.config_file):
        raise Exception(f"The primer3 config file {args.config_file} does not exist.")
    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)

    if args.min_overlap > args.max_overlap:
        raise Exception("The minimum overlap cannot be greater than the maximum overlap.")
    if args.min_overlap < 0:
        raise Exception("The minimum overlap cannot be less than 0.")
    if args.max_overlap < 0:
        raise Exception("The maximum overlap cannot be less than 0.")
    if args.overlap_increase < 0:
        raise Exception("The overlap increase cannot be less than 0.")
    if args.overlap_shift < 0:
        raise Exception("The overlap shift cannot be less than 0.")

    if args.amplicon_size < 0:
        raise Exception("The amplicon size cannot be less than 0.")
    if args.min_overlap > args.amplicon_size:
        raise Exception("The overlap cannot be greater than the amplicon size.")
    if args.max_overlap > args.amplicon_size:
        raise Exception("The overlap cannot be greater than the amplicon size.")
    
    # Load the regions
    regions = __load_regions(args.regions)

    # Load the sequence
    seq_record = __extract_sequence_record_from_fasta(args.fasta)

    # Generate the amplicons by each regions
    list_of_regions = [] 
    for i, row in regions.iterrows():
        if row["start"] > row["end"]:
            start = row["end"]
            end = row["start"]
        else:
            start = row["start"]
            end = row["end"]

        amplicons_with_primers = __generate_amplicon_and_primers_by_regions(
            region_name=row["loci"],
            region_start=start,
            region_end=end,
            amplicon_min_overlap=args.min_overlap,
            amplicon_max_overlap=args.max_overlap,
            amplicon_size=args.amplicon_size,
            overlap_increase_increment=args.overlap_increase,
            overlap_shift_increment=args.overlap_shift,
            primer3_settings_path=args.config_file,
            temp_dir=args.temp_dir,
            sequence=seq_record.seq
        )
        region = {
            "name": row["loci"],
            "amplicons": amplicons_with_primers
        }
        list_of_regions.append(region)
    
    # To Json output
    output = {
        "sequence": seq_record.id,
        "regions": list_of_regions
    }

    with open(args.output_file, "w") as f:
        json.dump(output, f, indent=4)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)