import sys
import asyncio
import logging
import subprocess
import re
import os

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

import pandas as pd

from collections import defaultdict
from configs import PrimerGenConfig
from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


class PrimerGenerator:
    """
    Generates the Primers for a given amplicon using Primer3
    """

    def __init__(
        self,
        region_name: str,
        amplicon_index: int,
        amplicon_sequence: Seq,
        pool_name: str,
        primer_ok_regions_list: tuple[int, int, int, int],
        config: PrimerGenConfig,
    ):
        self.amplicon_sequence = amplicon_sequence
        self.amplicon_id = f"{region_name}-{amplicon_index}"
        self.pool_name = pool_name
        self.primer_ok_regions_list = primer_ok_regions_list
        self._primer3_settings = config.primer3_settings
        self.temp_dir = config.temp_dir

    async def generate_primers(self) -> tuple[list | None, list | None]:
        """
        Generate the input file for primer3 and write it to a temporary file
        Take the output from primer 3 and parse it to this classes parsers
        Return the primer pairs along with additional information
        or None if either no forward or reverse primers
        """
        file_name = self.__write_temp_primer_gen_file()

        # Run primer3_core
        command = f"primer3_core < {file_name}"
        result = await asyncio.create_subprocess_shell(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        stdout, stderr = await result.communicate()

        if result.returncode != 0:
            raise Exception(f"Primer3 failed with error: {stderr}")

        return self.__parse_output_from_primer3(stdout)

    def __write_temp_primer_gen_file(self) -> str:
        PRIMER3_OK_REGIONS_LIST = f"{self.primer_ok_regions_list[0]},{self.primer_ok_regions_list[1]},{self.primer_ok_regions_list[2]},{self.primer_ok_regions_list[3]}"

        with open(os.path.join(self.temp_dir, self.amplicon_id), "w") as file:
            file.write(
                f"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST={PRIMER3_OK_REGIONS_LIST}\n"
            )
            file.write(f"SEQUENCE_ID={self.amplicon_id}\n")
            file.write(f"SEQUENCE_TEMPLATE={str(self.amplicon_sequence)}\n")
            file.write(self._primer3_settings)

            return file.name

    def __extract_primer_data(
        self, n_primers: int, output: str, left_or_right: str
    ) -> list | None:
        """Extract the primer data from the primer3 output. Extracts the sequence, tm, gc_percent and hairpin_th."""
        primers = []
        for i in range(0, n_primers):
            # Create a dictionary for each primer
            primer_entry = {}

            # Extract the index of the primer to later get the position
            raw_index_pattern = rf"PRIMER_{left_or_right}_{i}=(\d+),(\d+)\r?\n"
            index_pattern = re.compile(raw_index_pattern)
            index = re.search(index_pattern, output)
            primer_entry["index"] = int(index.group(1))

            # Extract the data for each primer
            raw_data_pattern = rf"PRIMER_{left_or_right}_{i}_(SEQUENCE|TM|GC_PERCENT|HAIRPIN_TH)=(\w+.*)\r?\n"
            data_pattern = re.compile(raw_data_pattern)
            data = re.findall(data_pattern, output)
            for entry in data:
                key = entry[0].lower()
                try:
                    value = float(entry[1])
                except ValueError:
                    value = entry[1]
                if key == "sequence":
                    key = "sequence"
                primer_entry[key] = value
            primers.append(primer_entry)

        if len(primers) == 0:
            return None

        return primers

    def __parse_output_from_primer3(
        self, output: bytes
    ) -> tuple[list | None, list | None]:
        """convert byte output to string output, afterwards find the desired primers and extract the data"""
        pattern = re.compile(r"PRIMER_(LEFT|RIGHT)_NUM_RETURNED=(\d+)\r?\n")

        output = str(output, "utf-8")
        pattern_search_result = re.findall(pattern, output)
        error_pattern = re.compile(rf"PRIMER_ERROR=[a-zA-Z0-9_\- ]+\r?\n")
        error_search_results = re.findall(error_pattern, output)
        if len(error_search_results) > 0:
            raise Exception(
                f"Primer3 failed with error(s):xw {','.join(error_search_results)}"
            )

        n_left_primers = int(pattern_search_result[0][1])
        n_right_primers = int(pattern_search_result[1][1])

        logging.info(
            f"Found {n_left_primers} left primers and {n_right_primers} right primers for sequence {self.amplicon_id} for pool {self.pool_name}."
        )

        # With the number of primers returned we iterate over the output
        # and extract the data for each primerpair. The data is stored in a dictionary
        # corresponding to the left and right primer. This is then stored in the overall results dictionary

        forward_primers = self.__extract_primer_data(n_left_primers, output, "LEFT")
        reverse_primers = self.__extract_primer_data(n_right_primers, output, "RIGHT")

        return forward_primers, reverse_primers


class RegionIterator:
    """Iterate over the regions in the database"""

    def __init__(self, regions: pd.DataFrame) -> None:
        self.regions = regions.iterrows()

    def __aiter__(self) -> "RegionIterator":
        return self

    async def __anext__(self) -> tuple[str, pd.Series]:
        try:
            region = next(self.regions)
        except StopIteration:
            raise StopAsyncIteration
        return region


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


def __clear_db(db: DBHandler) -> None:
    # clear the database proto_primers table
    db.execute("DELETE FROM proto_primers;")


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
    failed_amplicons: list[dict],
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
            failed_amplicons.append(
                {
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_name}",
                    "pool": pool_name,
                    "amplicon_start": adjusted_amplicon_start,
                    "amplicon_end": adjusted_amplicon_end,
                }
            )
            logging.warning(
                f"No primers generated for region {region_name} amplicon {idx}. Skipping amplicon."
            )
            continue

        # add the primers to the list
        for primer in forward_primers:
            list_of_primers.append(
                {
                    "pool": pool_name,
                    "region_name": region_name,
                    "amplicon_name": f"{region_name}-{idx}-{pool_name}",
                    "strand": "forward",
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
                    "strand": "reverse",
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

    # Clear the database
    __clear_db(db)

    # Setup parameters
    # pool_offset is the number of nucleotides that are skipped between each pool
    pool_offset = int((1 - config.min_overlap) * config.min_amplicon_size)
    # amplicon_offset is the number of nucleotides that are skipped between each amplicon in a pool
    amplicon_offset = int((1 - config.min_overlap * config.pool_count) * config.min_amplicon_size)
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
    failed_amplicons: list[dict] = []
    async for _, row in RegionIterator(regions=regions):
        start = row["start"]
        end = row["end"]
        name = row["name"]

        # Check if the region is valid
        if end > len(sequence_record.seq):
            logging.warning(
                f"Region {name} is out of bounds. Skipping region. Please check the region coordinates."
            )
            continue
        
        seeds = defaultdict(list[int])
        coords = defaultdict(list[tuple[int, int]])
        for pool in range(config.pool_count):
            # Generate the seed coordinates for each pool
            seeds[pool] = __generate_seed_coordinates(
                start + pool * pool_offset, end, config.min_amplicon_size, amplicon_offset
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
                failed_amplicons,
            )

    # Insert the primers into the database
    df = pd.DataFrame(list_of_primers)
    df.to_sql("proto_primers", db.conn, if_exists="append", index=False)

    # Insert the failed amplicons into the database
    df = pd.DataFrame(failed_amplicons)
    df.to_sql("failed_amplicons", db.conn, if_exists="append", index=False)


if __name__ == "__main__":
    # Run the main function using asyncio
    try:
        loop = asyncio.get_event_loop()
        loop.run_until_complete(main())
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
