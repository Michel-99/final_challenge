import eggnog_library as eggnog
import csv

# set up result file directory

result_1A = "results/1_A_homologs.txt"
result_1B = "results/1_B_unique_protein_IDs.txt"
result_1C = "results/1_C_functional_category_counts.csv"
result_detailed_2 = "results/2_detailed_results.txt"
result_3 = "results/3_universal_ogs.tsv"
summary_file = "results/COMPLETE_SUMMARY.txt"

# Dataframe setup
print("Setting up dataframes ...")
df_annotations = eggnog.dataframe_setup_annotations()
df_members = eggnog.dataframe_setup_members()
df_species = eggnog.dataframe_setup_taxid_info()
df_functional_categories_description = eggnog.dataframe_setup_functional_categories()

###----------------------------------------------------------------------
### 1) A) homologuous genes in humans and chimp but not mouse

# get species ids
human_id = eggnog.get_species_id_by_name("Homo sapiens", df_species)
chimp_id = eggnog.get_species_id_by_name("Pan troglodytes", df_species)
mouse_id = eggnog.get_species_id_by_name("Mus musculus", df_species)

# set which species should be included and excluded
include = [human_id, chimp_id]
exclude = [mouse_id]

# Scan for homologs within the inclusion list; omit excluded species. Assign return to variable "homologs"
print("Identify homologs in the target species set, filtering out excluded taxa ...")
homologs_df = eggnog.filter_by_ids(
    df_members, "species_taxid_containing_protein", include, exclude
)

# output to .txt file
print("Writing number of homologs to 1_A_homologs.txt file in /results ...")
with open(result_1A, "w") as f:
    f.write(
        f'We identified {len(homologs_df)} homologous genes shared by "{", ".join([eggnog.get_species_name_by_id(n, df_species) for n in include])}" that have diverged or are absent in "{", ".join([eggnog.get_species_name_by_id(n, df_species) for n in exclude])}"'
    )


###----------------------------------------------------------------------
### 1) B) extract unique protein IDs from homologs file and

print("Extracting unique protein IDs from homologs ...")

# extract unique protein IDs from homologs dataframe
unique_protein_ids = (
    homologs_df["protein_id"].str.split(",").explode().str.strip().unique()
)

print("Writing unique protein IDs to 1_B_unique_protein_IDs.txt file in /results ...")

with open(result_1B, "w") as f:
    f.write("\n".join(unique_protein_ids))

###----------------------------------------------------------------------
### 1) C) functional categories of homologous genes

# get functional categories for homologous genes by merging with annotations dataframe on orthologous group id
print("Extracting functional categories for homologous genes ...")

functional_categories_in_homologs = df_annotations[
    df_annotations["orthologous_group_id"].isin(
        homologs_df["orthologous_group_id"].tolist()
    )
]["functional_category"]

print("Counting occurrences of each functional category ...")

category_counts_df = (
    functional_categories_in_homologs.apply(list).explode().value_counts().reset_index()
)
category_counts_df.columns = ["category_code", "count"]

category_counts_df = category_counts_df.merge(
    df_functional_categories_description, on="category_code", how="left"
)

print("Exporting result table to results/1_C_functional_category_counts.csv ...")

category_counts_df.to_csv(result_1C, index=False)


###----------------------------------------------------------------------
### 1) D) ortholog genes found only in humans and chimp.


unique_orth_genes = homologs_df.loc[
    homologs_df["num_of_species"] == 2,
    ["orthologous_group_id", "species_taxid_containing_protein"],
]
print(f"{len(unique_orth_genes)} ortholog genes are only found in human and chimps")


###----------------------------------------------------------------------
### 1 E) ortholog genes found only in primates


primate_specific_OG = eggnog.filter_by_species_names(
    homologs_df, "species_taxid_containing_protein", eggnog.PRIMATES, df_species
)
print(f"{len(primate_specific_OG)} ortholog genes are only found in primates.")

###----------------------------------------------------------------------
### QUESTION 2: Lineage Analysis - Primates, Chicken, Fish vs Rodents
print("\n" + "=" * 80)
print("QUESTION 2: LINEAGE ANALYSIS")
print("Orthologs in Primates + Chicken + Fish, checking losses in Rodents")
print("=" * 80)

print("\nCleaning TaxIDs (removing protein suffixes)...")


# Apply the cleaning function
df_members["clean_taxid_set"] = df_members["species_taxid_containing_protein"].apply(
    eggnog.clean_taxid_string
)

print("Analyzing lineage conservation and loss...")

# Define Target IDs (verified these are standard TaxIDs)


# Retrieve sets
primates_q2 = eggnog.get_og_set(
    eggnog.IDS_EX2["human"], df_members
) | eggnog.get_og_set(eggnog.IDS_EX2["chimp"], df_members)
chicken = eggnog.get_og_set(eggnog.IDS_EX2["chicken"], df_members)
fish = eggnog.get_og_set(eggnog.IDS_EX2["danio"], df_members) | eggnog.get_og_set(
    eggnog.IDS_EX2["takifugu"], df_members
)
mouse = eggnog.get_og_set(eggnog.IDS_EX2["mouse"], df_members)
rat = eggnog.get_og_set(eggnog.IDS_EX2["rat"], df_members)

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

with open(result_detailed_2, "w") as f:
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

print(f"\n✅ Q2 Results saved to {result_detailed_2}")

###----------------------------------------------------------------------
### QUESTION 3: Universal Genes (99% or more of all animal species)

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

og_col = (
    "orthologous_group_id" if "orthologous_group_id" in universal_ogs.columns else None
)
if og_col:
    universal_ogs[[og_col, "actual_sp_count"]].to_csv(result_3, sep="\t", index=False)
else:
    universal_ogs[["actual_sp_count"]].to_csv(result_3, sep="\t")

print(f"✅ Q3 Results saved to {result_3}")

###----------------------------------------------------------------------
# save files

print("\n" + "=" * 80)
print("SAVING COMPREHENSIVE SUMMARY")
print("=" * 80)


with open(summary_file, "w") as f:
    f.write("=" * 80 + "\n")
    f.write("METAZOAN GENE CONSERVATION AND LOSS ANALYSIS - COMPLETE SUMMARY\n")
    f.write("Questions 1A-1E, Question 2, and Question 3\n")
    f.write("=" * 80 + "\n\n")

    f.write("QUESTION 1: PRIMATE-SPECIFIC GENE ANALYSIS\n")
    f.write("-" * 80 + "\n")
    f.write(f"1A) Homologs in Human & Chimp but NOT in Mouse: {len(homologs_df)} OGs\n")
    f.write(f"1B) Unique protein IDs: {len(unique_protein_ids)}\n")
    f.write(f"1C) Functional categories found: {len(category_counts_df)}\n")
    f.write(f"    Top 3 categories:\n")
    for idx, row in category_counts_df.head(3).iterrows():
        f.write(f"      - {row['category_code']}: {row['count']} genes\n")
    f.write(f"1D) OGs found ONLY in Human and Chimp: {len(unique_orth_genes)} OGs\n")
    f.write(f"1E) Primate-specific OGs: {len(primate_specific_OG)} OGs\n\n")

    f.write("QUESTION 2: LINEAGE ANALYSIS\n")
    f.write("-" * 80 + "\n")
    f.write(f"Core vertebrate OGs (Primates + Chicken + Fish): {len(core_set)}\n")
    f.write(f"Lost in BOTH Mouse and Rat: {len(lost_both)}\n")
    f.write(f"Lost ONLY in Mouse: {len(lost_only_mouse)}\n")
    f.write(f"Lost ONLY in Rat: {len(lost_only_rat)}\n\n")

    f.write("QUESTION 3: UNIVERSAL ANIMAL GENES\n")
    f.write("-" * 80 + "\n")
    f.write(f"Total unique species in dataset: {total_sp_count}\n")
    f.write(f"Universal OGs (99%+): {len(universal_ogs)}\n\n")

    f.write("=" * 80 + "\n")
    f.write("RESULT FILES SAVED:\n")
    f.write("-" * 80 + "\n")
    f.write("Question 1:\n")
    f.write("  - results/1_A_homologs.txt\n")
    f.write("  - results/1_B_unique_protein_IDs.txt\n")
    f.write("  - results/1_C_functional_category_counts.csv\n")
    f.write("Question 2:\n")
    f.write("  - results/2_detailed_results.txt\n")
    f.write("Question 3:\n")
    f.write("  - results/3_universal_ogs.tsv\n")
    f.write("=" * 80 + "\n")

print(f"\n✅ Complete summary saved to: {summary_file}")

print("\n" + "=" * 80)
print("ALL ANALYSES COMPLETE!")
print("=" * 80)
print(f"Main summary: {summary_file}")
print("All result files saved to: results/")
print("\n" + "=" * 80)
