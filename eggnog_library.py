"""Library for working with eggNOG files (eggNOG v5.0)"""

from pathlib import Path
from typing import List, Optional, Set
import pandas as pd
import re


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

IDS_EX2 = {
    "human": 9606,
    "chimp": 9598,
    "chicken": 9031,
    "danio": 7955,
    "takifugu": 31033,
    "mouse": 10090,
    "rat": 10116,
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

    Raises:
        FileNotFoundError: If the annotations file is not found in the data directory.
        pd.errors.ParserError: If the file is corrupted or malformed.
        pd.errors.EmptyDataError: If the file is empty.
    """
    try:
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
    except FileNotFoundError:
        raise FileNotFoundError(
            "File '33208_annotations.tsv' not found. Run: bash runall.sh to download files."
        )
    except pd.errors.ParserError as e:
        raise pd.errors.ParserError(
            f"File is corrupted or malformed. Re-download with: bash runall.sh\n"
            f"Details: {e}"
        ) from e
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError(
            "Data file is empty. Re-download with: bash runall.sh"
        ) from e
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

    Raises:
        FileNotFoundError: If the members file is not found in the data directory.
        pd.errors.ParserError: If the file is corrupted or malformed.
        pd.errors.EmptyDataError: If the file is empty.
    """
    try:
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
    except FileNotFoundError:
        raise FileNotFoundError(
            "File '33208_members.tsv' not found. Run: bash runall.sh to download files."
        )
    except pd.errors.ParserError as e:
        raise pd.errors.ParserError(
            f"File is corrupted or malformed. Re-download with: bash runall.sh\n"
            f"Details: {e}"
        ) from e
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError(
            "Data file is empty. Re-download with: bash runall.sh"
        ) from e
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

    Raises:
        FileNotFoundError: If the taxid info file is not found in the data directory.
        pd.errors.ParserError: If the file is corrupted or malformed.
        pd.errors.EmptyDataError: If the file is empty.
    """
    try:
        column_names = [
            "species_taxid",
            "species_name",
            "rank",
            "named_lineage",
            "tax_id_lineage",
        ]
        df = pd.read_csv(
            "data/e5.taxid_info.tsv", sep="\t", header=0, names=column_names
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            "File 'e5.taxid_info.tsv' not found. Run: bash runall.sh to download files."
        )
    except pd.errors.ParserError as e:
        raise pd.errors.ParserError(
            f"File is corrupted or malformed. Re-download with: bash runall.sh\n"
            f"Details: {e}"
        ) from e
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError(
            "Data file is empty. Re-download with: bash runall.sh"
        ) from e
    return df


def dataframe_setup_functional_categories() -> pd.DataFrame:
    """
    Set up the Dataframe for functional categories by parsing a text file.

    Returns:
        pd.DataFrame: A DataFrame containing functional categories with columns:
                      - Category_code
                      - Description

    Raises:
        FileNotFoundError: If the functional categories file is not found in the data directory.
        pd.errors.ParserError: If the file is corrupted or malformed.
        pd.errors.EmptyDataError: If the file is empty.
    """
    try:
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
    except FileNotFoundError:
        raise FileNotFoundError(
            "File 'eggnog4.functional_categories.txt' not found. Run: bash runall.sh to download files."
        )
    except pd.errors.ParserError as e:
        raise pd.errors.ParserError(
            f"File is corrupted or malformed. Re-download with: bash runall.sh\n"
            f"Details: {e}"
        )
    except pd.errors.EmptyDataError:
        raise pd.errors.EmptyDataError(
            "Data file is empty. Re-download with: bash runall.sh"
        )
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

    Raises:
        ValueError: If the species name is not found or has duplicate entries.

    Examples:
        >>> import pandas as pd
        >>> data = {'species_name': ['E. coli', 'H. sapiens'], 'species_taxid': [562, 9606]}
        >>> df = pd.DataFrame(data)
        >>> get_species_id_by_name('E. coli', df)
        562
    """
    try:
        species_id = df.loc[df["species_name"] == species_name, "species_taxid"].item()
    except ValueError:
        raise ValueError(
            f"Species '{species_name}' not found or has duplicate entries in dataframe."
        )
    return species_id


def get_species_name_by_id(species_id: int, df: pd.DataFrame) -> str:
    """
    Get the species name by species id from the species dataframe.

    Args:
        species_id (int): The species taxid.
        df (pd.DataFrame): The dataframe containing species information.

    Returns:
        str: The name of the species.

    Raises:
        ValueError: If the species ID is not found or has duplicate entries.

    Examples:
        >>> import pandas as pd
        >>> data = {'species_name': ['E. coli', 'H. sapiens'], 'species_taxid': [562, 9606]}
        >>> df = pd.DataFrame(data)
        >>> get_species_name_by_id(9606, df)
        'H. sapiens'
    """
    try:
        species_name = df.loc[df["species_taxid"] == species_id, "species_name"].item()
    except ValueError:
        raise ValueError(
            f"Species with ID '{species_id}' not found or has duplicate entries in dataframe."
        )
    return species_name


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

    Raises:
        KeyError: If the specified column is not found in the dataframe.

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

    try:
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
    except KeyError:
        raise KeyError(
            f"Column '{column}' not found in dataframe. "
            f"Available columns: {list(df.columns)}"
        )

    return df[mask]


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

    Raises:
        ValueError: If any species name in the iterable is not found in the species dataframe.
        KeyError: If the specified column is not found in the dataframe.

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

    try:
        species_ids = get_species_ids_from_names(species_names, df_species)
    except ValueError as e:
        raise ValueError(f"Failed to resolve species names to IDs: {e}") from e
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
    Raises:
        KeyError: If the specified column is not found in the dataframe.
    Example:
        >>> df = pd.DataFrame({'ids': ['1,2,3', '4,5', '2,3']})
        >>> allowed = ['1', '2', '3', '4', '5']
        >>> filter_allowed_ids(df, 'ids', allowed)
           ids
        0  1,2,3
        1  4,5
        2  2,3
    """
    allowed_set = set(str(item) for item in allowed_ids)

    try:
        mask = (
            df[column]
            .astype(str)
            .str.split(r"[,\s]+")
            .apply(lambda ids: set(ids) <= allowed_set)
        )
    except KeyError:
        raise KeyError(
            f"Column '{column}' not found in dataframe. "
            f"Available columns: {list(df.columns)}"
        )

    return df[mask]


def clean_taxid_string(raw_string: str) -> Set[int]:
    """
    Extract unique taxonomic IDs from a protein ID string.

    Parses a comma-separated string of protein IDs with taxonomic prefixes
    and returns a set of unique taxonomic IDs.

    Args:
        raw_string (str): Comma-separated protein IDs with taxid prefixes
                         (e.g., "9606.ENSP123, 10090.ENSMUSP456")

    Returns:
        Set[int]: Set of unique taxonomic IDs extracted from the string

    Raises:
        ValueError: If the string contains non-numeric taxonomic ID prefixes.

    Examples:
        >>> clean_taxid_string("9606.ENSP123, 10090.ENSMUSP456, 10116.ENSRNOP789")
        {9606, 10090, 10116}
        >>> clean_taxid_string("562.ECOLI_001")
        {562}
    """
    s = str(raw_string)
    items = s.split(",")
    try:
        clean_set = {int(item.split(".")[0].strip()) for item in items}
    except ValueError:
        raise ValueError(
            f"Failed to parse taxonomic IDs from string: '{raw_string}'. "
            f"Expected format: 'TAXID.PROTEINID' (e.g., '9606.ENSP00000123')."
        )
    return clean_set


def get_og_set(target_taxid: int, df_members: pd.DataFrame) -> Set[str]:
    """
    Get all orthologous group IDs containing a specific taxonomic ID.

    Filters the members dataframe to find all orthologous groups that contain
    the specified taxonomic ID. Requires that the dataframe has a 'clean_taxid_set'
    column created by applying clean_taxid_string.

    Args:
        target_taxid (int): The taxonomic ID to search for
        df_members (pd.DataFrame): Members dataframe with 'clean_taxid_set' column
                                   and 'orthologous_group_id' column or index

    Returns:
        Set[str]: Set of orthologous group IDs containing the target taxonomic ID

    Raises:
        KeyError: If the 'clean_taxid_set' column is missing from the dataframe.

    Examples:
        >>> import pandas as pd
        >>> df = pd.DataFrame({
        ...     'orthologous_group_id': ['OG1', 'OG2', 'OG3'],
        ...     'clean_taxid_set': [{9606, 10090}, {9606}, {10090, 10116}]
        ... })
        >>> get_og_set(9606, df)
        {'OG1', 'OG2'}
    """
    try:
        mask = df_members["clean_taxid_set"].apply(lambda s: target_taxid in s)
    except KeyError:
        raise KeyError(
            "Column 'clean_taxid_set' not found. Apply clean_taxid_string to the "
            "'species_taxid_containing_protein' column first."
        )
    og_col = (
        "orthologous_group_id"
        if "orthologous_group_id" in df_members.columns
        else df_members.index.name
    )
    if og_col in df_members.columns:
        return set(df_members.loc[mask, og_col])
    else:
        return set(df_members.loc[mask].index)
