import argparse
import pandas as pd
import json
import sys

DEFAULT_ADJACENCY_LIMIT = 500

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter primers by alignment")
    parser.add_argument(
        "--alignment",
        type=str,
        required=True,
        help="Alignment file with primer sequences",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file with primer filtered sequences",
    )
    parser.add_argument(
        "--primers",
        type=str,
        required=True,
        help="JSON file with original primer sequences",
    )
    parser.add_argument(
        "--adjacency_limit",
        type=int,
        default=DEFAULT_ADJACENCY_LIMIT,
        help=f"Limit within which primers of opposite strand are considered adjancent enough to be problematic. Default: {DEFAULT_ADJACENCY_LIMIT}",
    )

    return parser.parse_args()


def __write_output(primers: dict, output: str) -> None:
    with open(output, "w") as file:
        json.dump(primers, file, indent=4)


def __calculate_badness(alignment: pd.DataFrame, args: argparse.Namespace) -> list:
    """
    Takes in a dataframe containing the parsed output from bowtie
    and returns a list of primers to remove

    Algorithm:
        1. Find all primers not aligning to the original site:
            - Primers aligning multiple times
            - Primers with mismatches
            - Primers aligning to the wrong strand
        2. For each primer, find all primers aligning to the other strand within a certain distance not of the same amplicon
        3. Calculate a badness score for each primer considering:
            - the amount of mismatches
            - misalignments for the primer
            - the amount of adjacent primers (inverse the distance)
        4. Append the score in the filtered_proto_primers.json file and output
    """
    def calculate_badness(primer: pd.Series) -> float:
        """
        Calculates the badness score for a primer based on the amount of mismatches, misalignments and adjacent primers
        """
        badness = 1
        if primer["primer_strand"] == "forward":
            adjacency_limit = args.adjacency_limit
            opposite_strand = "reverse"
            adjacent_primers = alignment.loc[
                (alignment["primer_amplicon"] != primer["primer_amplicon"]) &
                (alignment["strand"] == opposite_strand) &
                (alignment["position"] <= primer["position"] + adjacency_limit) &
                (alignment["position"] >= primer["position"])
            ]
        if primer["primer_strand"] == "reverse":
            adjacency_limit = -1*args.adjacency_limit
            opposite_strand = "forward"
            adjacent_primers = alignment.loc[
                (alignment["primer_amplicon"] != primer["primer_amplicon"]) &
                (alignment["strand"] == opposite_strand) &
                (alignment["position"] >= primer["position"] + adjacency_limit) &
                (alignment["position"] <= primer["position"])
            ]

        print(f"Found {len(adjacent_primers)} adjacent primers for primer {primer['primer']}")
        badness += len(adjacent_primers)
        return badness
    
    problematic_primers = alignment.loc[
        (alignment["matches"] > 1) |
        (alignment["strand"] != alignment["primer_strand"]) |
        (alignment["mismatches_descriptor"].isna() == False)
    ].copy()
    problematic_primers.reset_index(inplace=True)

     # check each misaligned_primer for adjacent primers -> opposite strand, within adjacency_limit
    print(f"Found {len(problematic_primers)} misalignments")
    problematic_primers.loc[:,"badness"] = problematic_primers.apply(calculate_badness, axis=1)
    
    # join both tables
    alignment = alignment.merge(problematic_primers[["primer", "badness"]], on="primer", how="left", indicator=False)

    alignment.fillna(0, inplace=True)
    alignment = alignment.groupby(["primer_region", "primer_amplicon", "primer_strand", "primer_id"]).agg(
        badness=("badness", "sum"),
    ).reset_index()

    return alignment

def __add_values_to_json(original_primers: dict, alignment: pd.DataFrame, args: argparse.Namespace) -> dict:
    """
    Takes in a dataframe containing the parsed output from bowtie
    and a dictionary with the original primers and adds the badness score
    to each corresponding primer in the dictionary.
    """
    modified_primers = original_primers.copy()
    regions = original_primers["regions"]
    # For each region, amplicon and (forward|reverse) primer, add the badness score
    for i, region in enumerate(regions):
        region_name = region["region_name"]

        for j, amplicon in enumerate(region["amplicons"]):
            amplicon_name = amplicon["amplicon_name"]
            for k, primer in enumerate(amplicon["forward_primers"]):
                print(region_name, amplicon_name, "forward", k)
                primer["badness"] = alignment.loc[
                    (alignment["primer_region"] == region_name) &
                    (alignment["primer_amplicon"] == amplicon_name) &
                    (alignment["primer_strand"] == "forward") &
                    (alignment["primer_id"] == k)
                ]["badness"].values[0]
                amplicon["forward_primers"][k] = primer

            for k, primer in enumerate(amplicon["reverse_primers"]):
                primer["badness"] = alignment.loc[
                    (alignment["primer_region"] == region_name) &
                    (alignment["primer_amplicon"] == amplicon_name) &
                    (alignment["primer_strand"] == "reverse") &
                    (alignment["primer_id"] == k)
                ]["badness"].values[0]
                amplicon["reverse_primers"][k] = primer

            region["amplicons"][j] = amplicon
        modified_primers["regions"][i] = region
    
    return modified_primers

def main():
    args = get_args()

    # Read alignment file
    # Columns: primer	strand	chromosome	position	sequence	read_quality	matches	primer_region	primer_amplicon	primer_strand	primer_id
    dtypes = {
        "primer": "str",
        "strand": "str",
        "chromosome": "str",
        "position": "int",
        "sequence": "str",
        "read_quality": "str",
        "matches": "int",
        "primer_region": "str",
        "primer_amplicon": "str",
        "primer_strand": "str",
        "primer_id": "int",
    }
    alignment = pd.read_csv(args.alignment, sep="\t", header=0, dtype=dtypes)

    # Read primer file
    with open(args.primers, "r") as file:
        original_primers = json.load(file)

    # Filter primers by alignment

    # Write output
    alignment_with_badness = __calculate_badness(alignment, args)
    modified_primers = __add_values_to_json(original_primers, alignment_with_badness, args)
    __write_output(modified_primers, args.output)


if __name__ == "__main__":
    print("Starting primer filtering")
    main()
    try:
        pass
    except Exception as e:
        print(e)
        sys.exit(1)
