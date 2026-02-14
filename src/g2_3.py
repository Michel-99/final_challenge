import os
import sys
import pandas as pd

# Setup path
sys.path.append("/home/clara/Lectures/final_challenge")
os.chdir("/home/clara/Lectures/final_challenge")

import eggnog_library as eggnog
import csv

# Create results directory if it doesn't exist
if not os.path.exists("results"):
    os.makedirs("results")
if not os.path.exists("temp"):
    os.makedirs("temp")

print("=" * 80)
print("METAZOAN GENE CONSERVATION AND LOSS ANALYSIS")
print("Questions 1A-1E, Question 2, and Question 3")
print("=" * 80)

# Dataframe setup
print("\nSetting up dataframes ...")
df_annotations = eggnog.dataframe_setup_annotations()
df_members = eggnog.dataframe_setup_members()
df_species = eggnog.dataframe_setup_taxid_info()
df_functional_categories_description = eggnog.dataframe_setup_functional_categories()


###======================================================================
### QUESTION 2: Lineage Analysis - Primates, Chicken, Fish vs Rodents
###======================================================================

print("\n" + "=" * 80)
print("QUESTION 2: LINEAGE ANALYSIS")
print("Orthologs in Primates + Chicken + Fish, checking losses in Rodents")
print("=" * 80)

print("\nCleaning TaxIDs (removing protein suffixes)...")


def clean_taxid_string(raw_string):
    """
    Input: "9606.ENSP123, 10090.ENSMUSP456, 10116.ENSRNOP789"
    Output: Set{9606, 10090, 10116}
    """
    s = str(raw_string)
    items = s.split(",")
    clean_set = {int(item.split(".")[0].strip()) for item in items}
    return clean_set


# Apply the cleaning function
df_members["clean_taxid_set"] = df_members["species_taxid_containing_protein"].apply(
    clean_taxid_string
)

print("Analyzing lineage conservation and loss...")

# Define Target IDs (verified these are standard TaxIDs)
IDS = {
    "human": 9606,
    "chimp": 9598,
    "chicken": 9031,
    "danio": 7955,
    "takifugu": 31033,
    "mouse": 10090,
    "rat": 10116,
}


def get_og_set(target_taxid):
    # use mask instead string matching
    mask = df_members["clean_taxid_set"].apply(lambda s: target_taxid in s)
    og_col = (
        "orthologous_group_id"
        if "orthologous_group_id" in df_members.columns
        else df_members.index.name
    )
    if og_col in df_members.columns:
        return set(df_members.loc[mask, og_col])
    else:
        return set(df_members.loc[mask].index)


# Retrieve sets
primates_q2 = get_og_set(IDS["human"]) | get_og_set(IDS["chimp"])
chicken = get_og_set(IDS["chicken"])
fish = get_og_set(IDS["danio"]) | get_og_set(IDS["takifugu"])
mouse = get_og_set(IDS["mouse"])
rat = get_og_set(IDS["rat"])

print(f"\n=== SET SIZES ===")
print(f"Primates OGs: {len(primates_q2)}")
print(f"Chicken OGs: {len(chicken)}")
print(f"Fish OGs: {len(fish)}")
print(f"Mouse OGs: {len(mouse)}")
print(f"Rat OGs: {len(rat)}")

# Logic: Core Vertebrate Genes (present in Fish + Birds + Primates)
core_set = primates_q2 & chicken & fish

# Identify losses in Rodents
# (In Core Set) MINUS (Any Rodent)
lost_both = core_set - (mouse | rat)

# (In Core Set AND In Rat) MINUS (Mouse) -> Lost only in Mouse
lost_only_mouse = (core_set & rat) - mouse

# (In Core Set AND In Mouse) MINUS (Rat) -> Lost only in Rat
lost_only_rat = (core_set & mouse) - rat

print(f"\nCore vertebrate OGs (in Primates + Chicken + Fish): {len(core_set)}")
print(f"Lost in BOTH Mouse and Rat: {len(lost_both)}")
print(f"Lost ONLY in Mouse: {len(lost_only_mouse)}")
print(f"Lost ONLY in Rat: {len(lost_only_rat)}")

# Save Q2 Results
output_file_q2 = "results/2_detailed_results.txt"
with open(output_file_q2, "w") as f:
    f.write("Evolutionary Analysis: Primates, Chicken, Fish vs Rodents\n")
    f.write("=" * 60 + "\n")
    f.write(f"Conserved in Primates+Chicken+Fish: {len(core_set)}\n")
    f.write(f"Lost in both Mouse and Rat: {len(lost_both)}\n")
    f.write(f"Lost ONLY in Mouse: {len(lost_only_mouse)}\n")
    f.write(f"Lost ONLY in Rat: {len(lost_only_rat)}\n")
    f.write(
        "\nExample OGs lost in both:\n"
        + "\n".join(str(og) for og in list(lost_both)[:10])
    )

print(f"\n✅ Q2 Results saved to {output_file_q2}")

###======================================================================
### QUESTION 3: Universal Genes (99% or more of all animal species)
###======================================================================

print("\n" + "=" * 80)
print("QUESTION 3: UNIVERSAL ANIMAL GENES")
print("OGs present in 99% or more of all animal species")
print("=" * 80)

print("\nIdentifying universal animal genes...")

# 1. Calculate total species count from the data itself
all_species = set()
for s in df_members["clean_taxid_set"]:
    all_species.update(s)

total_sp_count = len(all_species)
threshold = total_sp_count * 0.99

print(f"Total unique species found in file: {total_sp_count}")
print(f"Threshold for universal (99%): {threshold:.2f}")

# 2. Count species per OG
df_members["actual_sp_count"] = df_members["clean_taxid_set"].apply(len)

# 3. Filter
universal_ogs = df_members[df_members["actual_sp_count"] >= threshold]

print(f"✅ Found {len(universal_ogs)} universal OGs (99%+ species)")

# Save Q3 Results
output_file_q3 = "results/3_universal_ogs.tsv"
og_col = (
    "orthologous_group_id" if "orthologous_group_id" in universal_ogs.columns else None
)
if og_col:
    universal_ogs[[og_col, "actual_sp_count"]].to_csv(
        output_file_q3, sep="\t", index=False
    )
else:
    universal_ogs[["actual_sp_count"]].to_csv(output_file_q3, sep="\t")

print(f"✅ Q3 Results saved to {output_file_q3}")
