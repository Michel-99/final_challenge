"""Library for working with eggNOG files (eggNOG v5.0)"""

from pathlib import Path
from typing import List, Optional, Set
import pandas as pd
import re
import csv

# --FUNCTIONS--
PRIMATES = {
    "Homo sapiens",
    "Tarsius syrichta",
    "Callithrix jacchus",
    "Macaca fascicularis",
    "Papio anubis",
    "Gorilla gorilla",
    "Pan paniscus",
    "Pan troglodytes",
    "Pongo abelii",
    "Saimiri boliviensis",
    "Chlorocebus sabaeus",
    "Rhinopithecus roxellana",
    "Nomascus leucogenys",
    "Otolemur garnettii",
    "Macaca mulatta",
}


# -- Dataframe setup functions--
def dataframe_setup_annotations() -> pd.DataFrame:
    """
    Set up the Dataframe for annotations. The orthologous_group_id is set as index.

    Returns:
        pd.DataFrame: A DataFrame containing annotation data with columns:
                      - evolutionary_level
                      - orthologous_group_id (index)
                      - functional category
                      - functional_description
    """
    column_names = [
        "evolutionary_level",
        "orthologous_group_id",
        "functional_category",
        "functional_description",
    ]
    df = pd.read_csv(
        "data/33208_annotations.tsv",
        sep="\t",
        header=None,
        names=column_names,
        # index_col="orthologous_group_id",
    )

    return df


def dataframe_setup_members() -> pd.DataFrame:
    """
    Set up the Dataframe for members. The orthologous_group_id is set as index.

    Returns:
        pd.DataFrame: A DataFrame containing member data with columns:
                      - evolutionary_level
                      - orthologous_group_id (index)
                      - num_of_proteins
                      - num_of_species
                      - protein_id
                      - species_taxid_containing_protein
    """
    column_names = [
        "evolutionary_level",
        "orthologous_group_id",
        "num_of_proteins",
        "num_of_species",
        "protein_id",
        "species_taxid_containing_protein",
    ]
    df = pd.read_csv(
        "data/33208_members.tsv",
        sep="\t",
        header=None,
        names=column_names,
        # index_col="orthologous_group_id",
    )

    return df


def dataframe_setup_taxid_info() -> pd.DataFrame:
    """
    Set up the Dataframe for tax (species) id info.

    Returns:
        pd.DataFrame: A DataFrame containing taxonomy info with columns:
                      - species_taxid
                      - species_name
                      - rank
                      - named_lineage
                      - tax_id_lineage
    """
    column_names = [
        "species_taxid",
        "species_name",
        "rank",
        "named_lineage",
        "tax_id_lineage",
    ]
    df = pd.read_csv("data/e5.taxid_info.tsv", sep="\t", header=0, names=column_names)
    return df


def dataframe_setup_functional_categories() -> pd.DataFrame:
    """
    Set up the Dataframe for functional categories by parsing a text file.

    Returns:
        pd.DataFrame: A DataFrame containing functional categories with columns:
                      - Category_code
                      - Description
    """
    data = []

    with open("data/eggnog4.functional_categories.txt", "r") as f:
        for line in f:
            line = line.strip()

            # Use regex to find lines starting with [Letter]
            # matches: [J] Description
            match = re.search(r"\[([A-Z])\]\s*(.*)", line)

            if match:
                letter = match.group(1)
                description = match.group(2)
                data.append([letter, description])

    df = pd.DataFrame(data, columns=["category_code", "description"])
    return df


# Functions using dataframes
def get_species_id_by_name(species_name: str, df: pd.DataFrame) -> int:
    """
    Get the species id by species name from the species dataframe.

    Args:
        species_name (str): The name of the species.
        df (pd.DataFrame): The dataframe containing species information.

    Returns:
        int: The species taxid.

    Examples:
        >>> import pandas as pd
        >>> data = {'species_name': ['E. coli', 'H. sapiens'], 'species_taxid': [562, 9606]}
        >>> df = pd.DataFrame(data)
        >>> get_species_id_by_name('E. coli', df)
        562
    """
    return df.loc[df["species_name"] == species_name, "species_taxid"].item()


def get_species_name_by_id(species_id: int, df: pd.DataFrame) -> str:
    """
    Get the species name by species id from the species dataframe.

    Args:
        species_id (int): The species taxid.
        df (pd.DataFrame): The dataframe containing species information.

    Returns:
        str: The name of the species.

    Examples:
        >>> import pandas as pd
        >>> data = {'species_name': ['E. coli', 'H. sapiens'], 'species_taxid': [562, 9606]}
        >>> df = pd.DataFrame(data)
        >>> get_species_name_by_id(9606, df)
        'H. sapiens'
    """
    return df.loc[df["species_taxid"] == species_id, "species_name"].item()


def get_species_ids_from_names(species_names, df_species: pd.DataFrame) -> List[int]:
    """
    Convert a list or set of species names to their corresponding IDs.

    Args:
        species_names: Iterable of species name strings
        df_species: Species dataframe containing taxid info

    Returns:
        List[int]: List of species IDs

    Examples:
        >>> import pandas as pd
        >>> data = {'species_name': ['E. coli', 'H. sapiens'], 'species_taxid': [562, 9606]}
        >>> df = pd.DataFrame(data)
        >>> get_species_ids_from_names(['E. coli', 'H. sapiens'], df)
        [562, 9606]
    """
    return [get_species_id_by_name(name, df_species) for name in species_names]


def filter_by_ids(
    df: pd.DataFrame,
    column: str,
    include_ids: List[str],
    exclude_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Filters a dataframe based on presence or absence of specific IDs
    within a comma-separated string column.

    Args:
        df (pd.DataFrame): The input dataframe.
        column (str): The column name containing comma-separated IDs.
        include_ids (List[str]): List of IDs that MUST be present.
        exclude_ids (Optional[List[str]]): List of IDs that MUST NOT be present. Defaults to None.

    Returns:
        pd.DataFrame: The filtered dataframe.

    Examples:
        >>> import pandas as pd
        >>> data = {'id': [1, 2, 3], 'tags': ['A,B', 'B,C', 'A,C']}
        >>> df = pd.DataFrame(data)
        >>> # Include 'A', Exclude 'B' -> Should match row 3 ('A,C')
        >>> filter_by_ids(df, 'tags', ['A'], ['B'])
           id tags
        2   3  A,C
        >>> # Include 'B' -> Should match row 1 and 2
        >>> filter_by_ids(df, 'tags', ['B'], [])
           id tags
        0   1  A,B
        1   2  B,C
    """
    if exclude_ids is None:
        exclude_ids = []

    # Start with a mask where everything is True
    mask = pd.Series(True, index=df.index)

    # 1. Must include ALL of these
    for inc_id in include_ids:
        # \b ensures we match "1" but not "11"
        regex_pattern = rf"\b{inc_id}\b"
        mask &= df[column].str.contains(regex_pattern, na=False)

    # 2. Must include NONE of these
    for exc_id in exclude_ids:
        regex_pattern = rf"\b{exc_id}\b"
        mask &= ~df[column].str.contains(regex_pattern, na=False)

    return df[mask]


# anpassen um pandas zu nutzen?


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


def filter_by_species_names(
    df: pd.DataFrame, column: str, species_names, df_species: pd.DataFrame
) -> pd.DataFrame:
    """
    Filter dataframe to include only rows where the specified column
    contains any of the given species (by name).

    This is a convenience wrapper that converts species names to IDs
    and then filters the dataframe.

    Args:
        df: Dataframe to filter
        column: Column name containing species IDs
        species_names: Iterable of species names to include
        df_species: Species dataframe for name-to-ID mapping

    Returns:
        pd.DataFrame: Filtered dataframe

    Examples:
        >>> import pandas as pd
        >>> df = pd.DataFrame({'ids': ['562,9606', '9606', '562']})
        >>> species_df = pd.DataFrame({
        ...     'species_name': ['E. coli', 'H. sapiens'],
        ...     'species_taxid': [562, 9606]
        ... })
        >>> filter_by_species_names(df, 'ids', ['E. coli', 'H. sapiens'], species_df)
           ids
        0  562,9606
        1  9606
        2  562
    """
    species_ids = get_species_ids_from_names(species_names, df_species)
    return filter_allowed_ids(df, column, species_ids)


def filter_allowed_ids(
    df: pd.DataFrame, column: str, allowed_ids: Set[str]
) -> pd.DataFrame:
    """
    Filter a DataFrame to include only rows where all IDs in a column are in the allowed set.
    This function splits the values in the specified column by common delimiters (commas and/or whitespace),
    then checks if all resulting IDs are present in the allowed_ids set. Only rows where all IDs are allowed
    are retained in the returned DataFrame.
    Args:
        df (pd.DataFrame): The input DataFrame to filter.
        column (str): The name of the column containing IDs to check. Values can be single IDs or
                      multiple IDs separated by commas and/or whitespace.
        allowed_ids (List[str]): A list of IDs that are considered allowed. Non-string IDs will be
                                 converted to strings for comparison.
    Returns:
        pd.DataFrame: A filtered DataFrame containing only rows where all IDs in the specified column
                      are present in the allowed_ids set. The original DataFrame is not modified.
    Example:
        >>> df = pd.DataFrame({'ids': ['1,2,3', '4,5', '2,3']})
        >>> allowed = ['1', '2', '3', '4', '5']
        >>> filter_allowed_ids(df, 'ids', allowed)
           ids
        0  1,2,3
        1  4,5
        2  2,3
    """
    allowed_set = set(str(id) for id in allowed_ids)

    mask = (
        df[column]
        .astype(str)
        .str.split(r"[,\s]+")
        .apply(lambda ids: set(ids) <= allowed_set)  # subset or equal
    )
    return df[mask]


# --END FUNCTIONS--
