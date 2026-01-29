import os
import csv


file_path = "/home/michel/Desktop/python/final_challenge/results/1_A_homologs_.tsv"
output_path = (
    "/home/michel/Desktop/python/final_challenge/results/1_B_unique_protein_IDs.txt"
)


def extract_protein_id(file_path: str) -> set:
    """
    Extract unique protein IDs from a tab-separated values (TSV) file.
    Reads a TSV file, skips the header row, and extracts protein IDs from the
    4th column (index 3). Handles comma-separated protein IDs within each cell.
    Args:
        file_path (str): The path to the TSV file to read.
    Returns:
        set: A set of unique protein IDs found in the file.
    Example:
        >>> protein_ids = extract_protein_id("data.tsv")
        >>> print(protein_ids)
        {'PROT001', 'PROT002', 'PROT003'}
    """

    unique_protein_ids = set()
    with open(file_path, "r") as f:
        rd = csv.reader(f, delimiter="\t")
        next(rd)  # Skip header row
        for line in rd:
            protein_ids = line[3].split(",")
            for protein_id in protein_ids:
                unique_protein_ids.add(protein_id)
    return unique_protein_ids


def write_protein_id(file_path: str, output_path: str) -> None:
    """
    Write unique protein IDs to a file, one per line.

    Args:
        output_path (str): The file path where protein IDs will be written.

    Returns:
        None
    """

    with open(output_path, "w") as file:
        unique = extract_protein_id(file_path)
        for id in unique:
            file.write(f"{id}\n")


write_protein_id(file_path, output_path)
