import os
import sys
import json
import asyncio

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd

from handler import PrimerGenerator, AmpliconGenerator
from config import Config


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
    pool_name: str,
    config: Config,
):
    """
    This function is called by __generate_amplicons_and_primers_in_pools and should
    not be called directly. It generates primers for a given amplicon using primer3.
    It it fails, then it returns (None, None) else it returns a predefined format
    in which the primers are present e.g.: (forward_primers, reverse_primers).
    """
    amplicon_forward_primers, amplicon_reverse_primers = [], []
    # the amplicon grows therefore we need to extract this and modify this value locally
    primer_ok_region_list = config.primer_ok_region_list
    # keep track of wether the amplicon is increasing in size too much into the pool_offset
    buffer_encroachment_counter  = 0
    while True:
        # Extract the amplicon sequence
        amplicon_sequence = sequence[start:end]
        # Try to generate primers for the amplicon
        primer_generator = PrimerGenerator(
            region_name=region_name,
            amplicon_index=idx,
            amplicon_sequence=amplicon_sequence,
            pool_name=pool_name,
            primer_ok_region_list=primer_ok_region_list,
            config=config,
        )
        forward_primers, reverse_primers = await primer_generator.generate_primers()

        amplicon_forward_primers.extend(forward_primers)
        amplicon_reverse_primers.extend(reverse_primers)

        if amplicon_forward_primers and amplicon_reverse_primers:
            break
        
        # prevent amplicons size to increase beyond 1/2 of the buffer region
        if buffer_encroachment_counter >= config.pool_offset / 2:
            break

        # Do not allow amplicons to be larger than max_amplicon_size
        if end - start > config.max_amplicon_size:
            break

        start -= config.amplicon_size_step
        end += config.amplicon_size_step
        primer_ok_region_list[0] += config.amplicon_size_step
        primer_ok_region_list[1] += config.amplicon_size_step
        buffer_encroachment_counter += config.amplicon_size_step

    return amplicon_forward_primers, amplicon_reverse_primers


def __split_region_into_pools(
    region_start: int,
    region_end: int,
    config: Config,
) -> tuple[list[dict], list[dict]]:
    """
    Divide a region into 2 pools and generate coordinates for potential amplicons.
    """
    pool_one_coords = []
    pool_two_coords = []
    amplicon_size = config.min_amplicon_size
    amplicon_offset = config.amplicon_offset
    pool_offset = config.pool_offset

    # add a buffer region to the start and end of the region, here: amplicon_size
    start = region_start - amplicon_size
    end = region_end + amplicon_size

    # pre_compute the start and ends of amplicons in the pool
    pool_one_pointer = start
    pool_two_pointer = start + pool_offset
    while pool_one_pointer < end:
        pool_one_coords.append(
            {
                "start": pool_one_pointer,
                "end": pool_one_pointer + amplicon_size,
            }
        )
        pool_one_pointer += amplicon_offset + amplicon_size

    while pool_two_pointer < end:
        pool_two_coords.append(
            {
                "start": pool_two_pointer,
                "end": pool_two_pointer + amplicon_size,
            }
        )
        pool_two_pointer += amplicon_offset + amplicon_size

    return pool_one_coords, pool_two_coords


async def __generate_primers_for_pool(
    region_name: str,
    pool_cords: list[dict],
    pool_name: str,
    sequence: Seq,
    config: Config,
):
    """Function to generate primers and extract final amplicons for a given pool and 'optimal' coordinates set."""

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
            pool_name=pool_name,
            config=config,
        )

        if not amplicon_forward_primers or not amplicon_reverse_primers:
            print(f"Failed to find primers for {region_name}-{idx} in region {coords['start']}-{coords['end']} in pool 1.")
            continue

        amplicon_forward_primers = __remove_duplicate_primers(amplicon_forward_primers)
        amplicon_reverse_primers = __remove_duplicate_primers(amplicon_reverse_primers)

        amplicons.append(
            {
                "amplicon_name": f"{region_name}-{idx}-{pool_name}",
                "forward_primers": amplicon_forward_primers,
                "reverse_primers": amplicon_reverse_primers,
            }
        )

    return amplicons


async def main():
    config = Config()

    # Load the regions
    regions = __load_regions(config.regions)

    # Load the sequence
    seq_record = __extract_sequence_record_from_fasta(config.fasta)

    # Generate pools for each region containing amplicons and primers
    pool_one_regions = []
    pool_two_regions = []
    async for i, row in AmpliconGenerator(regions=regions):
        if row["start"] > row["end"]:
            start = row["end"]
            end = row["start"]
        else:
            start = row["start"]
            end = row["end"]

        # Generate the optimal pool coordinates
        pool_one_coords, pool_two_coords = __split_region_into_pools(
            region_start=start, region_end=end, config=config
        )

        pool_one_amplicons_with_primers = await __generate_primers_for_pool(
            region_name=row["loci"],
            pool_cords=pool_one_coords,
            pool_name="pool_1",
            sequence=seq_record.seq,
            config=config,
        )
        pool_two_amplicons_with_primers = await __generate_primers_for_pool(
            region_name=row["loci"],
            pool_cords=pool_two_coords,
            pool_name="pool_2",
            sequence=seq_record.seq,
            config=config,
        )
        pool_one_regions.append(
            {"region_name": row["loci"], "amplicons": pool_one_amplicons_with_primers}
        )
        pool_two_regions.append(
            {"region_name": row["loci"], "amplicons": pool_two_amplicons_with_primers}
        )

    # To JSON
    pools = [
        {
            "pool_id": 1,
            "regions": pool_one_regions,
        },
        {
            "pool_id": 2,
            "regions": pool_two_regions,
        },
    ]
    with open(config.output_file, "w") as f:
        json.dump(pools, f, indent=4)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    # TODO: enable this
    # try:
    loop.run_until_complete(main())
    # except Exception as e:
    #     sys.stderr.write(f"ERROR: {e}\n")
    #     sys.exit(1)
