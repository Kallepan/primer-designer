
"""
Calculates proto primers for a given fasta file containing sequence. 
The primers are generated via primer3. Configuration options are defined in primer3_settings.conf.

"""
import os
import sys
import argparse
import yaml

import asyncio
import subprocess

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

    return parser

def _load_primer3_settings(config_file_path: str) -> dict:
    """
    Load the primer3 settings from the config file
    """
    settings = {}

    with open(config_file_path, "r") as handle:
        settings = yaml.safe_load(handle)
    
    settings_str = ""
    for k, v in settings.items():
        settings_str += f"{k}={v}\n"

    settings_str += "="

    return settings_str

async def _write_temp_seq_file(sequence_id: str, sequence_template: str, settings: str, amplicon_buffer_size: int) -> dict:
    """ Design the proto primers for the given sequence """
    # Calculate the product size range
    min_product_size = max(0, len(sequence_template) - amplicon_buffer_size)
        
    with open(f"tmp/{sequence_id}", "w") as temp_seq_file:
        temp_seq_file.write(f"PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{len(sequence_template)}")
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
    
    def __init__(self, primer3_settings: str, seq_file_path: str, amplicon_buffer_size: int):
        self._primer3_settings = primer3_settings
        self._seq_file = open(seq_file_path, "r")
        self._records = iter(SeqIO.parse(self._seq_file, "fasta"))
        self.amplicon_buffer_size = amplicon_buffer_size

    def __aiter__(self):
        return self

    def __del__(self):
        self._seq_file.close()

    async def __anext__(self):
        try:
            record = next(self._records)
        except StopIteration:
            raise StopAsyncIteration

        t = await _write_temp_seq_file(record.id , record.seq, self._primer3_settings, self.amplicon_buffer_size)
        
        return t
    
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
    
    settings_files = [path async for path in PrimerGenerator(primer3_settings, args.input, args.amplicon_buffer_size)]
    # TODO: iterate over each amplicon defined in the fasta file
    # implement multi threading to run primer3_core on each of the tmp files
    # take outpout and parse it into a dict
    # write the dict to a json file

if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    try:
        loop.run_until_complete(main())
    except StopAsyncIteration as e:
        sys.exit(0)
    except Exception as e:
        sys.stderr.write("ERROR: " + str(e))
        sys.exit(1)
