"""
Format the database to json to be used by the SADDLE SCRIPT
"""
import argparse
import pandas as pd
import logging
import json
import sys

from handlers import DBHandler

logging.basicConfig(level=logging.INFO)


def __get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format the database to json to be used by the SADDLE SCRIPT"
    )
    parser.add_argument(
        "--db", type=str, required=True, help="Path to the sqlite database"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to the output json file"
    )
    parser.add_argument(
        "--pool",
        type=str,
        required=True,
    )
    return parser.parse_args()


def __get_primer_object(primer: list) -> dict:
    """Returns a dictionary from the primer list taken from the database"""
    try:
        data = {
            "id": primer[0],
            "sequence": primer[5],
            "length": int(primer[6]),
            "tm": float(primer[7]),
            "gc_percent": float(primer[8]),
            "hairpin_th": float(primer[9]),
            "badness": float(primer[10]),
            "position": int(primer[11]),
        }
    except (TypeError, ValueError) as e:
        logging.error(f"Failed to parse primer {primer[0]}: {e}")
        # rescue if any of the values are null, should not happen
        data = {
            "id": primer[0],
            "sequence": primer[5],
            "length": primer[6],
            "tm": primer[7],
            "gc_percent": primer[8],
            "hairpin_th": primer[9],
            "badness": primer[10],
            "position": primer[11],
        }
    return data


def __get_primer_df(db: DBHandler, pool: str) -> pd.DataFrame:
    primers, column_names = db.select(
        """
            SELECT id, pool, region_name, amplicon_name, strand, sequence, length(sequence), tm, gc_percent, hairpin_th, badness, position
            FROM proto_primers
            WHERE 
                pool = ? AND 
                NOT (discarded)
            ORDER BY id ASC
        """,
        (pool,),
    )
    primer_df = pd.DataFrame(primers, columns=column_names, dtype="string")
    return primer_df


def __write_json(db: str, output: str, pool: str) -> None:
    primer_df = __get_primer_df(db, pool)

    # group dataframe by region then amplicon and finally then strand
    # Afterwards write to json
    json_dict = {"pool_id": pool, "regions": []}

    if primer_df.empty:
        with open(output, "w") as f:
            json.dump(json_dict, f, indent=4)
        return

    regions = primer_df.groupby("region_name")
    for region in regions:
        region_name = region[0]
        region_df = region[1]
        json_dict["regions"].append({"region_name": region_name, "amplicons": []})
        amplicons = region_df.groupby("amplicon_name")
        for amplicon in amplicons:
            amplicon_name = amplicon[0]
            amplicon_df = amplicon[1]
            json_dict["regions"][-1]["amplicons"].append(
                {
                    "amplicon_name": amplicon_name,
                    "forward_primers": [],
                    "reverse_primers": [],
                }
            )

            # Extract forward and reverse primers
            forward_primers = amplicon_df[amplicon_df["strand"] == "forward"]
            reverse_primers = amplicon_df[amplicon_df["strand"] == "reverse"]

            # Check if either forward or reverse primers are empty
            # If so, remove the amplicon
            if forward_primers.empty or reverse_primers.empty:
                json_dict["regions"][-1]["amplicons"].pop()
                continue

            for forward_primer in forward_primers.values:
                json_dict["regions"][-1]["amplicons"][-1]["forward_primers"].append(
                    __get_primer_object(forward_primer)
                )
            for reverse_primer in reverse_primers.values:
                json_dict["regions"][-1]["amplicons"][-1]["reverse_primers"].append(
                    __get_primer_object(reverse_primer)
                )

    with open(output, "w") as f:
        json.dump(json_dict, f, indent=4)


def main() -> None:
    logging.info("Formatting proto_primers table into json")
    args = __get_args()
    db = DBHandler(args.db)
    __write_json(db, args.output, args.pool)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
        sys.exit(1)
