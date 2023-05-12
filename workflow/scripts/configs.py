import os
import yaml
import argparse


def load_primer3_settings(path: str) -> str:
    """
    Static function used by the Config class.
    Load the primer3 settings from the config file.
    """
    KEYS_TO_IGNORE = [
        "PRIMER_PRODUCT_SIZE_RANGE",
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST",
    ]
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


class PrimerGenConfig:
    """
    Contains all the config for the proto primers script along with the arguments.
    """

    DEFAULT_MIN_AMPLICON_LENGTH = 200
    DEFAULT_MAX_AMPLICON_LENGTH = 300
    DEFAULT_AMPLICON_LENGTH_STEP = 5
    DEFAULT_AMPLICON_OFFSET = 50
    DEFAULT_POOL_OFFSET = round(DEFAULT_MIN_AMPLICON_LENGTH / 2)
    DEFAULT_FORWARD_END = 30
    DEFAULT_REVERSE_START = 30

    def __setup_args(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="""
            Generate a list of amplicons and generate Proto Primers for each amplicon from a sequence in a fasta file.
            Save the list of amplicons to a fasta file as well as the Proto Primers to a json.
            """
        )
        parser.add_argument(
            "-f",
            "--fasta",
            type=str,
            required=True,
            help="Path to the fasta file containing the sequence to be used to generate the amplicons.",
        )
        parser.add_argument(
            "-p",
            "--primer_3_settings_path",
            type=str,
            required=True,
            help="Path to the primer3 config file.",
        )
        parser.add_argument(
            "-r",
            "--regions",
            type=str,
            required=True,
            help="Path to the csv file containing the regions to be used to generate the amplicons.",
        )
        parser.add_argument(
            "-o",
            "--output_folder",
            type=str,
            required=True,
            help="Path to the output folder to be used to save the amplicons.",
        )
        parser.add_argument(
            "-t",
            "--temp_dir",
            type=str,
            help="Path to the temporary directory to be used when generating the amplicons.",
        )

        # Optional arguments
        parser.add_argument(
            "--pool_offset",
            type=int,
            default=self.DEFAULT_POOL_OFFSET,
            help=f"The amount by which the next pool is offset. Default is 1/2 default amplcon size {self.DEFAULT_POOL_OFFSET}.",
        )
        parser.add_argument(
            "--amplicon_offset",
            type=int,
            default=self.DEFAULT_AMPLICON_OFFSET,
            help=f"The amount by which the next amplicon in the same pool if offset. Default is {self.DEFAULT_AMPLICON_OFFSET}.",
        )
        parser.add_argument(
            "--min_amplicon_size",
            type=int,
            default=self.DEFAULT_MIN_AMPLICON_LENGTH,
            help=f"The minimum allowed size of the amplicons to be generated. Default is {self.DEFAULT_MIN_AMPLICON_LENGTH}.",
        )
        parser.add_argument(
            "--max_amplicon_size",
            type=int,
            default=self.DEFAULT_MAX_AMPLICON_LENGTH,
            help=f"The maximum allowed length of a generateed amplicon to be generated. Default is {self.DEFAULT_MAX_AMPLICON_LENGTH}.",
        )
        parser.add_argument(
            "--primer_ok_region_list",
            type=str,
            nargs=2,
            default=[self.DEFAULT_FORWARD_END, self.DEFAULT_REVERSE_START],
            help=f"Defines the bounds of the dynamically calculated area within the amplicon in which the forward and reverse primer must be located. The first index indicates the location of the forward primer (0 - <first_index>). The second index indicates the location of the reverse primer (len(amplicon)-<second_index> - len(amplicon)). Default is 0 to {self.DEFAULT_FORWARD_END} and len(amplicon)-{self.DEFAULT_REVERSE_START} to len(amplicon) respectively.",
        )
        parser.add_argument(
            "--amplicon_size_step",
            type=int,
            default=self.DEFAULT_AMPLICON_LENGTH_STEP,
            help=f"The amount by which the amplicon size is incremented on both ends. Default is {self.DEFAULT_AMPLICON_LENGTH_STEP}.",
        )
        return parser.parse_args()

    def __validate_args(self, args: argparse.Namespace) -> dict:
        # Error Handling
        if not os.path.exists(args.fasta):
            raise Exception(f"The fasta file {args.fasta} does not exist.")
        if not os.path.exists(args.regions):
            raise Exception(f"The regions file {args.regions} does not exist.")
        if not os.path.exists(args.primer_3_settings_path):
            raise Exception(
                f"The primer3 config file {args.primer_3_settings_path} does not exist."
            )
        if not os.path.exists(args.temp_dir):
            os.mkdir(args.temp_dir)

        # TODO: Add a check to make sure that the parameters are valid

        return vars(args)

    def __map_dict(self, args: dict) -> None:
        for k, v in args.items():
            setattr(self, k, v)

    def __init__(self):
        args = self.__setup_args()
        args_dict = self.__validate_args(args)
        self.__map_dict(args_dict)
        self.primer3_settings = load_primer3_settings(self.primer_3_settings_path)
