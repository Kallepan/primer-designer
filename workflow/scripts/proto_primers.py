import os
import sys
import json
import argparse
import yaml
import asyncio

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd

from handler import PrimerGenerator, AmpliconGenerator

DEFAULT_AMPLICON_SIZE=200
DEFAULT_OVERLAP_REGION_SIZE=50
DEFUALT_OVERLAP_REGION_INCREASE=10
DEFAULT_MAX_AMPLICON_LENGTH=250

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
        "--min_amplicon_size",
        type=int,
        default=DEFAULT_AMPLICON_SIZE,
        help=f"The minimum allowed size of the amplicons to be generated. Default is {DEFAULT_AMPLICON_SIZE}."
    )
    parser.add_argument(
        "--overlap_region_size",
        type=int,
        default=DEFAULT_OVERLAP_REGION_SIZE,
        help=f"The size of the overlap region to be used when generating the amplicons. Default is {DEFAULT_OVERLAP_REGION_SIZE}."
    )
    parser.add_argument(
        "--overlap_region_increase_increment",
        type=int,
        default=DEFUALT_OVERLAP_REGION_INCREASE,
        help=f"The increment to be used when increasing the overlap region size. Default is {DEFUALT_OVERLAP_REGION_INCREASE}."
    )
    parser.add_argument(
        "--max_amplicon_size",
        type=int,
        default=DEFAULT_MAX_AMPLICON_LENGTH,
        help=f"The maximum allowed length of a generateed amplicon to be generated. Default is {DEFAULT_MAX_AMPLICON_LENGTH}."
    )
    return parser

def __remove_duplicate_primers(list_of_primers: list[dict]) -> list:
    seen_primer = set()
    primers = []
    for primer_data in list_of_primers:
        if primer_data["primer_sequence"] in seen_primer:
            continue
        seen_primer.add(primer_data["primer_sequence"])
        primers.append(primer_data)
    return primers

async def __generate_primers(
    start: int,
    end: int,
    sequence: Seq,
    region_name: str,
    idx: int,
    primer3_settings: str,
    max_amplicon_size: int,
    overlap_region_increase_increment: int,
    overlap_region_size: int,
    pool_name: str,
    temp_dir: str,
    ):
    """
        This function is called by __generate_amplicons_and_primers_in_pools and should
        not be called directly. It generates primers for a given amplicon using primer3.
        It it fails, then it returns (None, None) else it returns a predefined format 
        in which the primers are present e.g.: (forward_primers, reverse_primers).
    """
    amplicon_forward_primers, amplicon_reverse_primers = [], []
    while True:
        amplicon_sequence = sequence[start:end]
        primer_generator = PrimerGenerator(
            region_name=region_name,
            amplicon_index=idx,
            amplicon_sequence=amplicon_sequence,
            temp_dir=temp_dir,
            primer3_settings=primer3_settings,
            overlap_region_size=overlap_region_size,
            pool_name=pool_name,
        )
        forward_primers, reverse_primers = await primer_generator.generate_primers()
        
        amplicon_forward_primers.extend(forward_primers)
        amplicon_reverse_primers.extend(reverse_primers)

        if amplicon_forward_primers and amplicon_reverse_primers:
            break;
        
        if end - start > max_amplicon_size:
            break;

        start -= overlap_region_increase_increment
        end += overlap_region_increase_increment
        overlap_region_size += overlap_region_increase_increment

    return amplicon_forward_primers, amplicon_reverse_primers

def __split_region_into_pools(
    region_start: int,
    region_end: int,
    min_amplicon_size: int,
    overlap_region_size: int,
) -> tuple[list[dict], list[dict]]:
    """Divide a region into 2 pools and generate coords for potential amplicons."""
    pool_one_coords = []
    pool_two_coords = []

    curr_pos = region_start - overlap_region_size
    end = region_end + overlap_region_size

    # pre_compute the start and ends of amplicons in the pool
    while curr_pos < end:
        pool_one_coords.append({
            "start": curr_pos,
            "end": curr_pos + min_amplicon_size, 
        })
        curr_pos += min_amplicon_size - overlap_region_size

        pool_two_coords.append({
            "start": curr_pos,
            "end": curr_pos + min_amplicon_size,
        })
        curr_pos += min_amplicon_size - overlap_region_size

    return pool_one_coords, pool_two_coords

async def __generate_primers_for_pool(
    region_name: str,
    pool_cords: list[dict],
    pool_name: str,
    sequence: Seq,
    max_amplicon_size: int,
    overlap_region_size: int,
    overlap_region_increase_increment: int,
    primer3_settings: str,
    temp_dir: str,
):
    """Function to generate primers and extract final amplicons for a given pool."""

    primer3_settings = __load_primer3_settings(primer3_settings)
    
    # generate the amplicons and primers for each pool
    # pool one
    amplicons = []
    for idx, coords in enumerate(pool_cords):
        amplicon_forward_primers, amplicon_reverse_primers = await __generate_primers(
            start=coords["start"],
            end=coords["end"],
            sequence=sequence,
            region_name=region_name,
            idx=idx,
            temp_dir=temp_dir,
            primer3_settings=primer3_settings,
            max_amplicon_size=max_amplicon_size,
            overlap_region_size=overlap_region_size,
            overlap_region_increase_increment=overlap_region_increase_increment,
            pool_name=pool_name,
        )

        if not amplicon_forward_primers or not amplicon_reverse_primers:
            print(f"Failed to find primers for {region_name}-{idx} in region {coords['start']}-{coords['end']} in pool 1.")
            continue
        
        amplicon_forward_primers = __remove_duplicate_primers(amplicon_forward_primers)
        amplicon_reverse_primers = __remove_duplicate_primers(amplicon_reverse_primers)

        amplicons.append({
            "amplicon_name": f"{region_name}-{idx}-{pool_name}",
            "forward_primers": amplicon_forward_primers,
            "reverse_primers": amplicon_reverse_primers,
        })

    return amplicons

async def main():
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

    if args.min_amplicon_size < 0:
        raise Exception("The amplicon size cannot be less than 0.")
    if args.overlap_region_size < 0:
        raise Exception("The buffer region cannot be less than 0.")
    if args.max_amplicon_size < 0:
        raise Exception("The max amplicon size cannot be less than 0.")

    # Load the regions
    regions = __load_regions(args.regions)

    # Load the sequence
    seq_record = __extract_sequence_record_from_fasta(args.fasta)

    # Generate pools for each region containing amplicons and primers
    pool_one = []
    pool_two = []
    async for i, row in AmpliconGenerator(regions=regions):
        if row["start"] > row["end"]:
            start = row["end"]
            end = row["start"]
        else:
            start = row["start"]
            end = row["end"]

        pool_one_coords, pool_two_coords = __split_region_into_pools(
            region_start=start,
            region_end=end,
            min_amplicon_size=args.min_amplicon_size,
            overlap_region_size=args.overlap_region_size,
        )

        pool_one_amplicons_with_primers = await __generate_primers_for_pool(
            region_name=row["loci"],
            pool_cords=pool_one_coords,
            pool_name="pool_1",
            sequence=seq_record.seq,
            max_amplicon_size=args.max_amplicon_size,
            overlap_region_size=args.overlap_region_size,
            overlap_region_increase_increment=args.overlap_region_increase_increment,
            primer3_settings=args.config_file,
            temp_dir=args.temp_dir,
        )
        pool_two_amplicons_with_primers = await __generate_primers_for_pool(
            region_name=row["loci"],
            pool_cords=pool_two_coords,
            pool_name="pool_2",
            sequence=seq_record.seq,
            max_amplicon_size=args.max_amplicon_size,
            overlap_region_size=args.overlap_region_size,
            overlap_region_increase_increment=args.overlap_region_increase_increment,
            primer3_settings=args.config_file,
            temp_dir=args.temp_dir,
        )
        pool_one.append({"region_name": row["loci"], "amplicons": pool_one_amplicons_with_primers})
        pool_two.append({"region_name": row["loci"], "amplicons": pool_two_amplicons_with_primers})
    
    # To JSON
    pools = {
        "pool_one": pool_one,
        "pool_two": pool_two,
    }
    with open(args.output_file, "w") as f:
        json.dump(pools, f, indent=4)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    # try:
    loop.run_until_complete(main())
    # except Exception as e:
    #     sys.stderr.write(f"ERROR: {e}\n")
    #     sys.exit(1)