import argparse
import sys
import sqlite3

def __parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--pool", type=int, required=True)
    return parser.parse_args()

def __setup_db_connection(input_db_file: str) -> sqlite3.Connection:
    con = sqlite3.connect(input_db_file)

    return con

def __write_fasta(db: sqlite3.Connection, args: argparse.Namespace) -> None:
    pools = db.execute(
        """
            SELECT DISTINCT pool FROM proto_primers ORDER BY pool ASC;
        """
    )

    with open(args.output, "w") as file:
        primers = db.execute(
            """
                SELECT primer_id, primer_sequence FROM proto_primers WHERE pool = ? ORDER BY primer_id ASC;
            """, 
            (args.pool,)
        ).fetchall()

        for primer in primers:
            id = primer[0]
            sequence = primer[1]
            file.write(f">{id}\n")
            file.write(sequence)
            file.write("\n")

def main():
    args = __parse_args()
    db = __setup_db_connection(args.input)
    __write_fasta(db, args)
    db.close()

if __name__ == "__main__":
    print("Formatting primers into fasta")
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)