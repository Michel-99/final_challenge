"""Library for working with eggNOG files (eggNOG v5.0)"""

from pathlib import Path
import pandas as pd
import re

# --FUNCTIONS--


# -- Dataframe setup functions--
def dataframe_setup_annotations():
    """Set up the Dataframe for annotations. The orthologous_group_id is set as index."""
    column_names = [
        "evolutionary_level",
        "orthologous_group_id",
        "functional category",
        "functional_description",
    ]
    df = pd.read_csv(
        "data/33208_annotations.tsv",
        sep="\t",
        header=None,
        names=column_names,
        index_col="orthologous_group_id",
    )

    return df


def dataframe_setup_members():
    """Set up the Dataframe for annotations. The orthologous_group_id is set as index."""
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
        index_col="orthologous_group_id",
    )

    return df


def dataframe_setup_taxid_info():
    """Set up the Dataframe for tax (species) id info. The species id is set as index."""
    column_names = [
        "species_taxid",
        "species_name",
        "rank",
        "named_lineage",
        "tax_id_lineage",
    ]
    df = pd.read_csv("data/e5.taxid_info.tsv", sep="\t", header=0, names=column_names)
    return df


def dataframe_setup_functional_categories():
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

    df = pd.DataFrame(data, columns=["Category_Letter", "Description"])
    return df


# Functions using dataframes
def get_species_id_by_name(species_name, df):
    """Get the species id by species name from the species dataframe."""
    species_id = df.loc[df["species_name"] == species_name, "species_taxid"].item()
    return species_id


def filter_by_ids(df, column, include_ids, exclude_ids):
    """
    Filters a dataframe based on presence or absence of specific IDs
    within a comma-separated string column.
    """
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
