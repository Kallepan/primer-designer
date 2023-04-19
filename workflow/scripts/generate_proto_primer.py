
"""
Calculates proto primers for a given fasta file containing sequence. 
The primers are generated via primer3. Configuration options are defined in primer3_settings.conf.

"""
import os
import sys
import argparse
import yaml

import asyncio
import re
import json

from Bio import SeqIO


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate primers for proto")
    parser.add_argument(
        "--input", "-i", help="Input fasta file to calculate the primers for", required=True, type=str
    )
    parser.add_argument(
        "--config", "-c", help="Path to the primer3 config file", required=True, type=str
    )
    
    parser.add_argument(
        "--amplicon_buffer_size", "-b", help="Buffer size for amplicons", required=True, type=int
    )

    parser.add_argument(
        "--temp_dir", "-t", help="Path to the temp directory", required=True, type=str
    )

    parser.add_argument(
        "--output_file", "-o", help="Path to the output file", required=True, type=str
    )

    return parser

def _load_primer3_settings(config_file_path: str) -> dict:
    """
    Load the primer3 settings from the config file
    """
    KEYS_TO_IGNORE = ["PRIMER_PRODUCT_SIZE_RANGE", "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"]
    settings = {}

    with open(config_file_path, "r") as handle:
        settings = yaml.safe_load(handle)
    
    settings_str = ""

    for k in KEYS_TO_IGNORE:
        settings.pop(k, None)

    for k, v in settings.items():
        settings_str += f"{k}={v}\n"

    settings_str += "="

    return settings_str

async def _write_temp_seq_file(sequence_id: str, sequence_template: str, settings: str, amplicon_buffer_size: int, temp_dir: str) -> dict:
    """ Design the proto primers for the given sequence """
    # Calculate the product size range
    PRIMER_3_OK_REGION_LIST=f"0,{amplicon_buffer_size},{len(sequence_template) - amplicon_buffer_size},{amplicon_buffer_size}"
    with open(os.path.join(temp_dir, str(sequence_id)), "w") as temp_seq_file:
        temp_seq_file.write(f"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST={PRIMER_3_OK_REGION_LIST}\n")
        temp_seq_file.write(f"SEQUENCE_ID={sequence_id}\n")
        temp_seq_file.write(f"SEQUENCE_TEMPLATE={str(sequence_template)}\n")
        temp_seq_file.write(settings)
        return temp_seq_file.name

class PrimerGenerator():
    """
    Generates the Primers in an async manner
    It needs the settings for primer3 as well as the file path to the fasta file
    and the buffer size for the amplicons
    """

    def __init__(self, primer3_settings: str, seq_file_path: str, amplicon_buffer_size: int, tmp_dir: str):
        self._primer3_settings = primer3_settings
        self._seq_file = open(seq_file_path, "r")
        self._records = iter(SeqIO.parse(self._seq_file, "fasta"))
        self._amplicon_buffer_size = amplicon_buffer_size
        self.tmp_dir = tmp_dir

    def __aiter__(self):
        return self

    def __del__(self):
        self._seq_file.close()

    async def __anext__(self):
        try:
            record = next(self._records)
        except StopIteration:
            raise StopAsyncIteration
        
        file_name = await _write_temp_seq_file(record.id , record.seq, self._primer3_settings, self._amplicon_buffer_size, self.tmp_dir)
        
        # Run primer3_core
        command = f"primer3_core < {file_name}"
        result = await asyncio.create_subprocess_shell(
            command,
            shell=True,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await result.communicate()
        
        if result.returncode != 0:
            raise Exception(f"Primer3 failed with error: {stderr}")

        return str(record.id), self.__parse_output_from_primer3(stdout, record.id) 
    
    def __extract_primer_data(self, n_primers: int, output: str, left_or_right: str) -> list:
        primers = []
        for i in range(1, n_primers):
            pattern = rf"PRIMER_{left_or_right}_{i}_(SEQUENCE|TM|GC_PERCENT|HAIRPIN_TH)=(\w+.*)\r?\n"
            data_pattern = re.compile(pattern)
            data = re.findall(data_pattern, output)
            primer_data = {}
            for entry in data:
                primer_data[entry[0]] = entry[1]
                primer_data["ID"] = f"{i}_{left_or_right}"
            primers.append(primer_data)
        return primers
    
    def __parse_output_from_primer3(self, output: bytes, sequence_id: str) -> dict:
        # convert byte output to string output, afterwards find the desired sequences
        pattern = re.compile(r"PRIMER_(LEFT|RIGHT)_NUM_RETURNED=(\d+)\r?\n")
        result = {}
        output = str(output, "utf-8")
        pattern_search_result = re.findall(pattern, output)

        if pattern_search_result is None or len(pattern_search_result) != 2:
            raise Exception("Primer3 failed to return any primers")
        
        n_left_primers = int(pattern_search_result[0][1])
        n_right_primers = int(pattern_search_result[1][1])

        print(f"Found {n_left_primers} left primers and {n_right_primers} right primers for sequence {sequence_id}")

        """
        With the number of primers returned we iterate over the output
        and extract the data for each primerpair. The data is stored in a dictionary
        corresponding to the left and right primer. This is then stored in the overall results dictionary
        """

        result["forward_primers"] = self.__extract_primer_data(n_left_primers, output, "LEFT")
        result["reverse_primers"] = self.__extract_primer_data(n_right_primers, output, "RIGHT")

        return result

async def main():
    parser = get_parser()
    args = parser.parse_args()

    # Read in the fasta files
    if not os.path.exists(args.input):
        raise Exception(f"Input file {args.input} does not exist")
    
    # Read in the primer3 settings
    if not os.path.exists(args.config):
        raise Exception("Primer3 config yaml does not exist")

    primer3_settings = _load_primer3_settings(args.config)
    
    # The magic happens here
    amplicons = []
    async for id, result in PrimerGenerator(primer3_settings, args.input, args.amplicon_buffer_size, args.temp_dir):
        primers = {
            "id": id,
            "forward_primers": result["forward_primers"],
            "reverse_primers": result["reverse_primers"]
        }
        amplicons.append(primers)
    
    with open(args.output_file, "w") as output_file:
        json.dump(amplicons, output_file, indent=4)

if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    try:
        loop.run_until_complete(main())
    except StopAsyncIteration as e:
        sys.exit(0)
    except Exception as e:
        sys.stderr.write("ERROR: " + str(e))
        sys.exit(1)
