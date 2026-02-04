import os
import sys
import pandas as pd
from collections import Counter

import sys

sys.path.append("/home/clara/Lectures/final_challenge")
import eggnog_library as eggnog


def run_analysis():
    print(f"Running from: {os.getcwd()}")
    print("Loading data into DataFrames...")

    # Load members data
    df_members = eggnog.dataframe_setup_members()

    # Standardize 6 columns (handles the index issue)
    df_members = df_members.reset_index()
    df_members.columns = [
        "taxid_level",
        "orthologous_group",
        "num_proteins",
        "num_species",
        "protein_ids",
        "species_taxid_containing_protein",
    ]

    # IMPORTANT: Ensure TaxID column is treated as a list of strings, not a single string
    # This prevents counting characters instead of species
    df_members["taxid_list"] = (
        df_members["species_taxid_containing_protein"]
        .astype(str)
        .apply(lambda x: x.split(","))
    )

    # question 2
    print("Analyzing lineage conservation and loss")

    # TaxIDs for target groups
    IDS = {
        "human": "9606",
        "chimp": "9598",
        "chicken": "9031",
        "danio": "7955",
        "takifugu": "31033",
        "mouse": "10090",
        "rat": "10116",
    }

    def get_og_set(taxid):
        # Check if the taxid exists in our newly created list of IDs
        mask = df_members["taxid_list"].apply(lambda x: taxid in x)
        return set(df_members.loc[mask, "orthologous_group"])

    # Build the Sets
    primates = get_og_set(IDS["human"]) | get_og_set(IDS["chimp"])
    chicken = get_og_set(IDS["chicken"])
    fish = get_og_set(IDS["danio"]) | get_og_set(IDS["takifugu"])
    mouse = get_og_set(IDS["mouse"])
    rat = get_og_set(IDS["rat"])

    # Intersect: Shared by Primates AND Chicken AND Fish
    core_set = primates & chicken & fish

    # Rodent Losses
    lost_both = core_set - (mouse | rat)
    lost_only_mouse = (core_set & rat) - mouse
    lost_only_rat = (core_set & mouse) - rat

    # Save Results
    with open("results/question_2_detailed_results.txt", "w") as f:
        f.write("Evolutionary Analysis: Primates, Chicken, Fish vs Rodents\n")
        f.write("=" * 60 + "\n")
        f.write(f"Conserved in Primates+Chicken+Fish: {len(core_set)}\n")
        f.write(f"Lost in both Mouse and Rat: {len(lost_both)}\n")
        f.write(f"Lost ONLY in Mouse: {len(lost_only_mouse)}\n")
        f.write(f"Lost ONLY in Rat: {len(lost_only_rat)}\n")
        f.write("\nExample OGs lost in both:\n" + "\n".join(list(lost_both)[:10]))

    # Question 3

    print("Identifying universal animal genes (q3)...")

    # Count total unique species mentioned in the entire file
    all_species = set()
    for row in df_members["taxid_list"]:
        all_species.update(row)

    total_sp_count = len(all_species)
    threshold = total_sp_count * 0.99
    print(f"Total species: {total_sp_count}, 99% Threshold: {threshold:.2f}")

    # Calculate count for each OG and filter
    df_members["actual_sp_count"] = df_members["taxid_list"].apply(len)
    universal_ogs = df_members[df_members["actual_sp_count"] >= threshold]

    # Save q3
    universal_ogs[["orthologous_group", "actual_sp_count"]].to_csv(
        "results/q3_universal_ogs.tsv", sep="\t", index=False
    )

    print(f"Analysis complete! Universal OGs found: {len(universal_ogs)}")


if __name__ == "__main__":
    # Ensure results directory exists
    if not os.path.exists("results"):
        os.makedirs("results")
    run_analysis()
