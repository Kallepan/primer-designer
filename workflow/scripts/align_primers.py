import argparse
import subprocess
import sys
import logging

DEFAULT_NUMBER_OF_MISMATCHES = 3

logging.basicConfig(level=logging.INFO)


def __get_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align primers to reference genome")

    parser.add_argument(
        "--primers", type=str, required=True, help="Fasta file containing primers"
    )
    parser.add_argument("--index", type=str, required=True, help="Path to index files")
    parser.add_argument(
        "--output", type=str, required=False, help="Path to output file"
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=DEFAULT_NUMBER_OF_MISMATCHES,
        help=f"Number of mismatches allowed. Default: {DEFAULT_NUMBER_OF_MISMATCHES}",
    )

    return parser.parse_args()


def __run_bowtie(args: argparse.Namespace):
    shell_cmd = f"bowtie -v {args.mismatches} -a -x {args.index} -f {args.primers}"
    sp = subprocess.run(
        shell_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
    )
    if sp.returncode != 0:
        raise Exception(sp.stderr.decode("utf-8"))

    return sp.stdout.decode("utf-8")


def main():
    logging.info("Aligning primers to a given genome")

    args = __get_parser()
    raw_alignment = __run_bowtie(args)

    with open(args.output, "w") as f:
        f.write(raw_alignment)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e, file=sys.stderr)
        sys.exit(1)
