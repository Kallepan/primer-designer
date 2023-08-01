import argparse
import sys
import logging

from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


def __parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--pool", type=int, required=True)
    return parser.parse_args()


def __write_fasta(db: DBHandler, args: argparse.Namespace) -> None:
    with open(args.output, "w") as file:
        primers, _ = db.select(
            """
                SELECT id, sequence FROM proto_primers WHERE pool = ? AND discarded = 0 ORDER BY id ASC;
            """,
            (args.pool,),
        )
        for primer in primers:
            id = primer[0]
            sequence = primer[1]
            file.write(f">{id}\n")
            file.write(sequence)
            file.write("\n")


def main():
    logging.info("Formatting primers into fasta")
    args = __parse_args()
    db = DBHandler(args.db)
    __write_fasta(db, args)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
