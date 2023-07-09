import os
import yaml
import argparse


class PrimerGenConfig:
    """
    Contains all the config for the proto primers script along with the arguments.
    """

    DEFAULT_MIN_AMPLICON_SIZE = 200
    DEFAULT_MAX_AMPLICON_SIZE = 300
    DEFAULT_MIN_OVERLAP = 0.2

    def __setup_args(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="""
            Generate a list of amplicons and generate Proto Primers for each amplicon from a sequence in a fasta file.
            Save the list of amplicons to a fasta file as well as the Proto Primers to a json.
            """
        )
        parser.add_argument(
            "-f",
            "--fasta_path",
            type=str,
            required=True,
            help="Path to the fasta file containing the sequence to be used to generate the amplicons.",
        )
        parser.add_argument(
            "-d",
            "--db_path",
            type=str,
            required=True,
            help="Path to the sqlite database",
        )
        parser.add_argument(
            "-p",
            "--primer_3_settings_path",
            type=str,
            required=True,
            help="Path to the primer3 config file.",
        )
        parser.add_argument(
            "-t",
            "--temp_dir",
            required=True,
            type=str,
            help="Path to the temporary directory to be used when generating the amplicons.",
        )

        # Optional arguments
        parser.add_argument(
            "--min_overlap",
            type=float,
            default=self.DEFAULT_MIN_OVERLAP,
            help=f"The minimum amount of overlap of two amplicons between pools in decimal. Default is {self.DEFAULT_MIN_OVERLAP}.",
        )
        parser.add_argument(
            "--min_amplicon_size",
            type=int,
            default=self.DEFAULT_MIN_AMPLICON_SIZE,
            help=f"The minimum allowed size of the amplicons to be generated. Default is {self.DEFAULT_MIN_AMPLICON_SIZE}.",
        )
        parser.add_argument(
            "--max_amplicon_size",
            type=int,
            default=self.DEFAULT_MAX_AMPLICON_SIZE,
            help=f"The maximum allowed length of a generateed amplicon to be generated. Default is {self.DEFAULT_MAX_AMPLICON_SIZE}.",
        )

        return parser.parse_args()

    def __validate_args(self, args: argparse.Namespace) -> dict:
        # Validate the arguments and return a dict of the arguments
        # Error Handling
        if not os.path.exists(args.fasta_path):
            raise Exception(f"The fasta file {args.fasta_path} does not exist.")
        if not os.path.exists(args.primer_3_settings_path):
            raise Exception(
                f"The primer3 config file {args.primer_3_settings_path} does not exist."
            )
        if not os.path.exists(args.temp_dir):
            os.mkdir(args.temp_dir)

        return vars(args)

    def __map_dict(self, args: dict) -> None:
        # For each key in the dict, set the value as an attribute of the class
        for k, v in args.items():
            setattr(self, k, v)

    @staticmethod
    def load_primer3_settings(path: str) -> str:
        """Static function used by the Config class. Load the primer3 settings from the config file."""
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

    def __init__(self):
        args = self.__setup_args()
        args_dict = self.__validate_args(args)
        self.__map_dict(args_dict)
        self.primer3_settings = self.load_primer3_settings(self.primer_3_settings_path)
