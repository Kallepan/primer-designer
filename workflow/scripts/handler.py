import re
import subprocess
import os
import logging
import asyncio

from Bio.Seq import Seq
from configs import PrimerGenConfig


class PrimerGenerator:
    """
    Generates the Primers for a given amplicon using Primer3
    """

    def __init__(
        self,
        region_name: str,
        amplicon_index: int,
        amplicon_sequence: Seq,
        pool_id: int,
        primer_ok_region_list: tuple[int, int],
        config: PrimerGenConfig,
    ):
        self.amplicon_sequence = amplicon_sequence
        self.amplicon_id = f"{region_name}-{amplicon_index}"
        self.pool_id = pool_id
        self.primer_ok_region_list = primer_ok_region_list
        self._primer3_settings = config.primer3_settings
        self.temp_dir = config.temp_dir

    async def generate_primers(self) -> tuple[list, list]:
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
        PRIMER3_OK_REGION_LIST = f"0,{self.primer_ok_region_list[0]},{len(self.amplicon_sequence)-self.primer_ok_region_list[1]},{self.primer_ok_region_list[1]}"
        with open(os.path.join(self.temp_dir, self.amplicon_id), "w") as file:
            file.write(
                f"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST={PRIMER3_OK_REGION_LIST}\n"
            )
            file.write(f"SEQUENCE_ID={self.amplicon_id}\n")
            file.write(f"SEQUENCE_TEMPLATE={str(self.amplicon_sequence)}\n")
            file.write(self._primer3_settings)

            return file.name

    def __extract_primer_data(
        self, n_primers: int, output: str, left_or_right: str
    ) -> list:
        primers = []
        for i in range(0, n_primers):
            pattern = rf"PRIMER_{left_or_right}_{i}_(SEQUENCE|TM|GC_PERCENT|HAIRPIN_TH)=(\w+.*)\r?\n"
            data_pattern = re.compile(pattern)
            data = re.findall(data_pattern, output)
            primer_entry = {}
            for entry in data:
                key = entry[0].lower()
                try:
                    value = float(entry[1])
                except ValueError:
                    value = entry[1]
                if key == "sequence":
                    key = "sequence"
                primer_entry[key] = value
            primer_entry["length"] = len(primer_entry["sequence"])
            primers.append(primer_entry)
        return primers

    def __parse_output_from_primer3(self, output: bytes) -> tuple[list, list]:
        # convert byte output to string output, afterwards find the desired sequences
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
            f"Found {n_left_primers} left primers and {n_right_primers} right primers for sequence {self.amplicon_id} for pool {self.pool_id}."
        )

        """
        With the number of primers returned we iterate over the output
        and extract the data for each primerpair. The data is stored in a dictionary
        corresponding to the left and right primer. This is then stored in the overall results dictionary
        """

        forward_primers = self.__extract_primer_data(n_left_primers, output, "LEFT")
        reverse_primers = self.__extract_primer_data(n_right_primers, output, "RIGHT")

        return forward_primers, reverse_primers


import pandas as pd
import asyncio


class AmpliconGenerator:
    def __init__(self, regions: pd.DataFrame):
        self.regions = regions.iterrows()

    def __aiter__(self):
        return self

    async def __anext__(self) -> tuple[str, pd.Series]:
        try:
            region = next(self.regions)
        except StopIteration:
            raise StopAsyncIteration
        return region
