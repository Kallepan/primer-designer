import os
import sys
import asyncio
import logging

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd
import numpy as np

from handler import PrimerGenerator, AmpliconGenerator
from configs import PrimerGenConfig
from db import DBHandler

logging.basicConfig(level=logging.INFO)


def __load_regions_from_db(db: DBHandler) -> pd.DataFrame:
    df = pd.read_sql_query(
        "SELECT * FROM regions",
        db.conn,
        dtype={
            "name": "string",
            "start": "int64",
            "end": "int64",
            "sequence": "string",
        },
    )
    return df


def __load_regions(path_to_file: str) -> pd.DataFrame:
    df = pd.read_json(
        path_to_file,
        dtype={
            "name": "string",
            "start": "int64",
            "end": "int64",
        },
    )

    # switch up the start and end columns if the start is greater than the end
    df["start"], df["end"] = np.where(
        df["start"] > df["end"], (df["end"], df["start"]), (df["start"], df["end"])
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
        if primer_data["sequence"] in seen_primer:
            continue
        seen_primer.add(primer_data["sequence"])
        primers.append(primer_data)
    return primers


async def __generate_primers(
    start: int,
    end: int,
    sequence: Seq,
    region_name: str,
    idx: int,
    pool_id: int,
    config: PrimerGenConfig,
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
    buffer_encroachment_counter = 0
    while True:
        # Extract the amplicon sequence
        amplicon_sequence = sequence[start:end]
        # Try to generate primers for the amplicon
        primer_generator = PrimerGenerator(
            region_name=region_name,
            amplicon_index=idx,
            amplicon_sequence=amplicon_sequence,
            pool_id=pool_id,
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
    config: PrimerGenConfig,
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
    pool_id: str,
    sequence: Seq,
    config: PrimerGenConfig,
    list_of_primers: list[dict],
):
    """Function to generate primers and extract final amplicons for a given pool and 'optimal' coordinates set."""

    # generate the amplicons and primers for each pool
    # keep track of those in pd.Dataframe object to later export to sqlite database
    for idx, coords in enumerate(pool_cords):
        amplicon_forward_primers, amplicon_reverse_primers = await __generate_primers(
            start=coords["start"],
            end=coords["end"],
            sequence=sequence,
            region_name=region_name,
            idx=idx,
            pool_id=pool_id,
            config=config,
        )

        if not amplicon_forward_primers or not amplicon_reverse_primers:
            logging.warn(
                f"Failed to find primers for {region_name}-{idx} in region {coords['start']}-{coords['end']} in {pool_id}."
            )
            continue

        amplicon_forward_primers = __remove_duplicate_primers(amplicon_forward_primers)
        amplicon_reverse_primers = __remove_duplicate_primers(amplicon_reverse_primers)

        for primer in amplicon_forward_primers:
            list_of_primers.append(
                {
                    "pool": pool_id,
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_id}",
                    "strand": "forward",
                    "sequence": primer["sequence"],
                    "length": primer["length"],
                    "tm": primer["tm"],
                    "gc_percent": primer["gc_percent"],
                    "hairpin_th": primer["hairpin_th"],
                    "badness": 0.0,
                }
            )

        for primer in amplicon_reverse_primers:
            list_of_primers.append(
                {
                    "pool": pool_id,
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_id}",
                    "strand": "reverse",
                    "sequence": primer["sequence"],
                    "length": primer["length"],
                    "tm": primer["tm"],
                    "gc_percent": primer["gc_percent"],
                    "hairpin_th": primer["hairpin_th"],
                    "badness": 0.0,
                }
            )

    return list_of_primers


async def main():
    config = PrimerGenConfig()
    # Load the regions
    db = DBHandler(config.db)
    regions = __load_regions_from_db(db=db)

    # Load the sequence
    seq_record = __extract_sequence_record_from_fasta(config.fasta)

    # Generate pools for each region containing amplicons and primers
    list_of_primers = []
    async for i, row in AmpliconGenerator(regions=regions):
        start = row["start"]
        end = row["end"]

        if start > end:
            raise ValueError(f"Start {start} is greater than end {end}.")
        if start < 0:
            raise ValueError(f"Start {start} is less than 0.")
        if end > len(seq_record.seq):
            raise ValueError(
                f"End {end} is greater than sequence length {len(seq_record.seq)}."
            )

        # Generate the optimal pool coordinates
        pool_one_coords, pool_two_coords = __split_region_into_pools(
            region_start=start, region_end=end, config=config
        )

        await __generate_primers_for_pool(
            region_name=row["name"],
            pool_cords=pool_one_coords,
            pool_id=0,
            sequence=seq_record.seq,
            config=config,
            list_of_primers=list_of_primers,
        )

        await __generate_primers_for_pool(
            region_name=row["name"],
            pool_cords=pool_two_coords,
            pool_id=1,
            sequence=seq_record.seq,
            config=config,
            list_of_primers=list_of_primers,
        )

    # Export the primers to a sqlite database
    df = pd.DataFrame(list_of_primers)
    df.to_sql("proto_primers", db.conn, if_exists="append", index=False)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
    try:
        pass
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
