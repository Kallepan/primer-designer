"""
Format the database to json to be used by the SADDLE SCRIPT
"""
import argparse
import pandas as pd
import json
from db import DBHandler

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Format the database to json to be used by the SADDLE SCRIPT")
    parser.add_argument(
        "--db",
        type=str,
        required=True,
        help="Path to the sqlite database"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output json file"
    )
    parser.add_argument(
        "--pool",
        type=str,
        required=True,
    )
    return parser.parse_args()

def __write_json(db: str, output: str, pool: str) -> None:
    primers = db.select(
        """
            SELECT id, pool, region_name, amplicon_name, strand, sequence, length, tm, gc_percent, hairpin_th, badness
            FROM proto_primers
            WHERE pool = ?
            ORDER BY id ASC
        """, (pool,))
    
    df = pd.DataFrame(
        primers, 
        columns=["id", "pool", "region_name", "amplicon_name", "strand", "sequence", "length", "tm", "gc_percent", "hairpin_th", "badness"],
        dtype="string"
    )
    
    # group dataframe by region then amplicon and finally then strand
    # Afterwards write to json
    json_dict = {
        "pool_id": pool,
        "regions": [

        ]
    }

    if df.empty:
        with open(output, "w") as f:
            json.dump(json_dict, f, indent=4)
        return

    regions = df.groupby("region_name")
    for region in regions:
        region_name = region[0]
        region_df = region[1]
        json_dict["regions"].append({
            "region_name": region_name,
            "amplicons": []
        })
        amplicons = region_df.groupby("amplicon_name")
        for amplicon in amplicons:
            amplicon_name = amplicon[0]
            amplicon_df = amplicon[1]
            json_dict["regions"][-1]["amplicons"].append({
                "amplicon_name": amplicon_name,
                "forward_primers": [],
                "reverse_primers": []
            })
            forward_primers = amplicon_df[amplicon_df["strand"] == "forward"]
            reverse_primers = amplicon_df[amplicon_df["strand"] == "reverse"]
            for forward_primer in forward_primers.values:
                json_dict["regions"][-1]["amplicons"][-1]["forward_primers"].append({
                    "id": forward_primer[0],
                    "sequence": forward_primer[5],
                    "length": forward_primer[6],
                    "tm": forward_primer[7],
                    "gc_percent": forward_primer[8],
                    "hairpin_th": forward_primer[9],
                    "badness": forward_primer[10]
                })
            for reverse_primer in reverse_primers.values:
                json_dict["regions"][-1]["amplicons"][-1]["reverse_primers"].append({
                    "id": reverse_primer[0],
                    "sequence": reverse_primer[5],
                    "length": reverse_primer[6],
                    "tm": reverse_primer[7],
                    "gc_percent": reverse_primer[8],
                    "hairpin_th": reverse_primer[9],
                    "badness": reverse_primer[10]
                })

    with open(output, "w") as f:
        json.dump(json_dict, f, indent=4)

def main() -> None:
    print("Formatting proto_primers table into json")
    args = get_args()
    db = DBHandler(args.db)
    __write_json(db, args.output, args.pool)

if __name__ == "__main__":
    main()
