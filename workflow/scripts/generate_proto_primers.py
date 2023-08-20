import asyncio
import logging

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd

from primer_gen.configs import PrimerGenConfig
from primer_gen.generator import PrimerGenerator
from primer_gen.utils import RegionIterator
from collections import defaultdict
from handlers import DBHandler

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


def __extract_sequence_record_from_fasta(path: str) -> SeqRecord:
    with open(path, "r") as file:
        seq = list(SeqIO.parse(file, "fasta"))

    if len(seq) == 0:
        raise Exception("No sequences found in fasta file")

    if len(seq) > 1:
        raise Exception("More than one sequence found in fasta file")

    return seq[0]


def __generate_seed_coordinates(
    start: int, end: int, amplicon_size: int, amplicon_offset: int
) -> list[int]:
    """Generates the seed coordinates for a given region. Each seed coordinate leads to an amplicon."""
    seed_coordinates = range(start, end, amplicon_offset + amplicon_size)
    return seed_coordinates


def __generate_amplicon_coordinates(
    seed_coordinates: list[int], amplicon_size: int
) -> list[tuple[int, int]]:
    """Generates the amplicon's start and end coordinates for a given list of seed coordinates."""
    amplicon_coordinates = []
    for seed in seed_coordinates:
        start = int(seed - amplicon_size / 2)
        end = int(seed + amplicon_size / 2)
        amplicon_coordinates.append((start, end))

    return amplicon_coordinates


async def __generate_primers(
    region_name: str,
    pool_name: str,
    amplicon_coords: list[tuple[int, int]],
    sequence: Seq,
    primer_ok_regions_list: list[int, int, int, int],
    config: PrimerGenConfig,
    list_of_primers: list[dict],
    amplicons: list[dict],
) -> None:
    """Generate the primer pairs for a given list of amplicon coordinates."""

    for idx, (amplicon_start, amplicon_end) in enumerate(amplicon_coords):
        # Ensure valid amplicon coordinates
        adjusted_amplicon_start = max(0, amplicon_start)
        adjusted_amplicon_end = min(len(sequence), amplicon_end)

        # extract the amplicon sequence
        amplicon_sequence = sequence[adjusted_amplicon_start:adjusted_amplicon_end]

        # create the primer generator
        primer_generator = PrimerGenerator(
            region_name=region_name,
            amplicon_index=idx,
            amplicon_sequence=amplicon_sequence,
            pool_name=pool_name,
            primer_ok_regions_list=primer_ok_regions_list,
            config=config,
        )

        # generate the primers
        forward_primers, reverse_primers = await primer_generator.generate_primers()

        # check if forward AND reverse primers were generated.
        if forward_primers is None or reverse_primers is None:
            amplicons.append(
                {
                    "region_name": region_name,
                    "name": f"{region_name}-{idx}-{pool_name}",
                    "pool": pool_name,
                    "start": adjusted_amplicon_start,
                    "end": adjusted_amplicon_end,
                    "failed": True,
                }
            )
            logging.warning(
                f"No primers generated for region {region_name} amplicon {idx}. Skipping amplicon."
            )
            continue

        # add the amplicon to the list
        amplicons.append(
            {
                "region_name": region_name,
                "name": f"{region_name}-{idx}-{pool_name}",
                "pool": pool_name,
                "start": adjusted_amplicon_start,
                "end": adjusted_amplicon_end,
                "failed": False,
            }
        )

        # add the primers to the list
        for primer in forward_primers:
            list_of_primers.append(
                {
                    "pool": pool_name,
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_name}",
                    "strand": "+",
                    "sequence": primer["sequence"],
                    "tm": primer["tm"],
                    "gc_percent": primer["gc_percent"],
                    "hairpin_th": primer["hairpin_th"],
                    "position": primer["index"] + adjusted_amplicon_start,
                }
            )

        for primer in reverse_primers:
            list_of_primers.append(
                {
                    "pool": pool_name,
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_name}",
                    "strand": "-",
                    "sequence": primer["sequence"],
                    "tm": primer["tm"],
                    "gc_percent": primer["gc_percent"],
                    "hairpin_th": primer["hairpin_th"],
                    "position": primer["index"] + adjusted_amplicon_start,
                }
            )


async def main():
    # Due to plenty of options, we use a config class to store the parameters
    config = PrimerGenConfig()

    # Load the regions
    db = DBHandler(config.db_path)
    regions = __load_regions_from_db(db)

    # Load the sequence
    sequence_record = __extract_sequence_record_from_fasta(config.fasta_path)

    # Setup parameters
    # pool_offset is the number of nucleotides that are skipped between each pool
    pool_offset = int((1 - config.min_overlap) * config.min_amplicon_size)
    # amplicon_offset is the number of nucleotides that are skipped between each amplicon in a pool
    amplicon_offset = int(
        (1 - config.min_overlap * config.pool_count) * config.min_amplicon_size
    )
    # if pool_count is 1, set amplion_offset to 0
    amplicon_offset = 0 if config.pool_count == 1 else amplicon_offset

    # amplicon_buffer is the number of nucleotides that are added to the amplicon size to allow for primer placement
    amplicon_buffer = int((config.max_amplicon_size - config.min_amplicon_size) / 2)
    # Primer ok regions are regions where the primer can be placed without exceeding the amplicon size
    # forward primers are placed in the range from 0 to amplicon_buffer
    # reverse primers are placed in the ragen from max_amplicon_size - amplicon_buffer to max_amplicon_size

    primer_ok_regions_list = [
        0,
        amplicon_buffer,
        config.max_amplicon_size - amplicon_buffer,
        amplicon_buffer,
    ]

    # Create the RegionIterator
    list_of_primers: list[dict] = []
    amplicons: list[dict] = []
    async for _, row in RegionIterator(regions=regions):
        start = row["start"]
        end = row["end"]
        name = row["name"]

        # Check if the region is valid
        if start < 0:
            logging.warning(
                f"Region {name} is out of bounds. Skipping region. Please check the region coordinates."
            )
            continue
        if end > len(sequence_record.seq):
            logging.warning(
                f"Region {name} is out of bounds. Skipping region. Please check the region coordinates."
            )
            continue
        if end - start < config.min_amplicon_size:
            logging.warning(
                f"Region {name} is too small. Skipping region. Please check the region coordinates."
            )
            continue

        seeds = defaultdict(list[int])
        coords = defaultdict(list[tuple[int, int]])
        for pool in range(config.pool_count):
            # Generate the seed coordinates for each pool
            seeds[pool] = __generate_seed_coordinates(
                start + pool * pool_offset,
                end,
                config.min_amplicon_size,
                amplicon_offset,
            )

            # Generate the start and end coordinates for each amplicon
            coords[pool] = __generate_amplicon_coordinates(
                seeds[pool], config.max_amplicon_size
            )

            # Generate the primers for each amplicon
            await __generate_primers(
                name,
                str(pool),
                coords[pool],
                sequence_record.seq,
                primer_ok_regions_list,
                config,
                list_of_primers,
                amplicons,
            )

    # Insert the primers into the database
    df = pd.DataFrame(list_of_primers)
    df.to_sql("proto_primers", db.conn, if_exists="append", index=False, chunksize=1000)

    # Insert the amplicons into the database
    df = pd.DataFrame(amplicons)
    df.to_sql("amplicons", db.conn, if_exists="append", index=False, chunksize=1000)


if __name__ == "__main__":
    # Run the main function using asyncio
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
