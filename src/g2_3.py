import os
import sys
import pandas as pd

# 1. Setup path
sys.path.append("/home/clara/Lectures/final_challenge")
import eggnog_library as eggnog


def run_analysis():
    print(f"Running from: {os.getcwd()}")
    print("Loading data into DataFrames...")

    # Load members data using your custom library
    df_members = eggnog.dataframe_setup_members()

    # Reset index to ensure we have a flat dataframe
    df_members = df_members.reset_index()
    df_members.columns = [
        "taxid_level",
        "orthologous_group",
        "num_proteins",
        "num_species",
        "protein_ids",
        "species_taxid_containing_protein",
    ]

    # change --> ID cleaning, did not find IDs before
    print("Cleaning TaxIDs (removing protein suffixes)...")

    def clean_taxid_string(raw_string):
        """
        Input: "9606.ENSP123, 10090.ENSMUSP456, 10116.ENSRNOP789"
        Output: Set{'9606', '10090', '10116'}
        """
        # 1. Convert to string to handle any NaN values safely
        s = str(raw_string)

        # 2. Split by comma to get individual items
        items = s.split(",")

        # 3. Clean each item:
        #    - split('.')[0] takes the number before the dot
        #    - strip() removes any accidental spaces
        clean_set = {item.split(".")[0].strip() for item in items}
        return clean_set

    # Apply the cleaning function
    df_members["clean_taxid_set"] = df_members[
        "species_taxid_containing_protein"
    ].apply(clean_taxid_string)

    # --- QUESTION 2: Lineage Analysis ---
    print("Analyzing lineage conservation and loss...")

    # Define Target IDs (verified these are standard TaxIDs)
    IDS = {
        "human": "9606",
        "chimp": "9598",
        "chicken": "9031",
        "danio": "7955",
        "takifugu": "31033",
        "mouse": "10090",
        "rat": "10116",
    }

    def get_og_set(target_taxid):
     #use mask instead string matching --> https://stackoverflow.com/questions/32280556/how-to-filter-a-dataframe-column-of-lists-for-those-that-contain-a-certain-item
        mask = df_members["clean_taxid_set"].apply(lambda s: target_taxid in s)
        return set(df_members.loc[mask, "orthologous_group"])

    # Retrieve sets
    primates = get_og_set(IDS["human"]) | get_og_set(IDS["chimp"])
    chicken = get_og_set(IDS["chicken"])
    fish = get_og_set(IDS["danio"]) | get_og_set(IDS["takifugu"])
    mouse = get_og_set(IDS["mouse"])
    rat = get_og_set(IDS["rat"])

    # Logic: Core Vertebrate Genes (present in Fish + Birds + Primates)
    core_set = primates & chicken & fish

    # Identify losses in Rodents
    # (In Core Set) MINUS (Any Rodent)
    lost_both = core_set - (mouse | rat)

    # (In Core Set AND In Rat) MINUS (Mouse) -> Lost only in Mouse
    lost_only_mouse = (core_set & rat) - mouse

    # (In Core Set AND In Mouse) MINUS (Rat) -> Lost only in Rat
    lost_only_rat = (core_set & mouse) - rat

    # Save Q2 Results
    output_file_q2 = "results/question_2_detailed_results.txt"
    with open(output_file_q2, "w") as f:
        f.write("Evolutionary Analysis: Primates, Chicken, Fish vs Rodents\n")
        f.write("=" * 60 + "\n")
        f.write(f"Conserved in Primates+Chicken+Fish: {len(core_set)}\n")
        f.write(f"Lost in both Mouse and Rat: {len(lost_both)}\n")
        f.write(f"Lost ONLY in Mouse: {len(lost_only_mouse)}\n")
        f.write(f"Lost ONLY in Rat: {len(lost_only_rat)}\n")
        f.write("\nExample OGs lost in both:\n" + "\n".join(list(lost_both)[:10]))

    print(
        f"Q2 Done. Found {len(core_set)} core genes. Results saved to {output_file_q2}"
    )

    # --- QUESTION 3: Universal Genes ---
    print("Identifying universal animal genes (q3)...")

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

    # Save Q3 Results
    output_file_q3 = "results/q3_universal_ogs.tsv"
    universal_ogs[["orthologous_group", "actual_sp_count"]].to_csv(
        output_file_q3, sep="\t", index=False
    )

    print(f"Q3 Analysis complete! Universal OGs found: {len(universal_ogs)}")


if __name__ == "__main__":
    if not os.path.exists("results"):
        os.makedirs("results")
    run_analysis()
