from Bio import SeqIO
import pandas as pd

def extract_seq_record_from_fasta_file(path_to_fasta_file: str) -> SeqIO.SeqRecord:
    """Read a fasta file and return a list of sequences"""

    with open(path_to_fasta_file, "r") as handle:
        fa = list(SeqIO.parse(handle, "fasta"))
    
    if len(fa) == 0:
        raise Exception("No sequences found in fasta file")
    
    if len(fa) > 1:
        raise Exception("More than one sequence found in fasta file")

    return fa[0]

def load_loci(path_to_loci_file: str) -> pd.DataFrame:
    """Read a loci csv file and return a list of loci"""

    df = pd.read_csv(
        path_to_loci_file, 
        sep=",", 
        header=0, 
        names=["loci", "start", "end"], 
        dtype={"loci": str, "start": int, "end": int}
    )

    return df

def extract_region(seq, start: int, end: int, buffer: int = 50) -> SeqIO.SeqRecord:
    """Extract a region of genome from a sequence"""
    
    if start > end:
        return seq[end-buffer:start+buffer]

    return seq[start-buffer:end+buffer]

def split_region_into_amplicons(region, amplicon_size: int = 1000, amplicon_buffer: int = 50) -> list[SeqIO.SeqRecord]:
    """Split a region of the genome into amplicons"""

    amplicons = []
    for i in range(0, len(region), amplicon_size):
        amplicons.append(region[i:i+amplicon_size])

    return amplicons

def write_amplicons_to_file(amplicons : dict, path_to_output_file: str) -> None:
    """Write amplicons to file"""

    with open(path_to_output_file, "w") as handle:
        for region_name, region_amplicons in amplicons.items():
            for i, amplicon in enumerate(region_amplicons):
                handle.write(f">{region_name}_{i}")
                handle.write("\n")
                handle.write(f"{amplicon}")
                handle.write("\n")