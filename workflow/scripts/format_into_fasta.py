import argparse
import sys

from db import DBHandler

def __parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--pool", type=int, required=True)
    return parser.parse_args()

def __write_fasta(db: DBHandler, args: argparse.Namespace) -> None:
    pools = db.execute(
        """
            SELECT DISTINCT pool FROM proto_primers ORDER BY pool ASC;
        """
    )

    with open(args.output, "w") as file:
        primers, column_names = db.select(
            """
                SELECT id, sequence FROM proto_primers WHERE pool = ? ORDER BY id ASC;
            """,
            (args.pool,)
        )
        for primer in primers:
            id = primer[0]
            sequence = primer[1]
            file.write(f">{id}\n")
            file.write(sequence)
            file.write("\n")

def main():
    print("Formatting primers into fasta")
    args = __parse_args()
    db = DBHandler(args.db)
    __write_fasta(db, args)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)