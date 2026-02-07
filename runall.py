import eggnog_library as eggnog
import csv

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
with open("results/1_A_homologs.txt", "w") as f:
    f.write(
        f'We identified {len(homologs_df)} homologous genes shared by "{", ".join([eggnog.get_species_name_by_id(n, df_species) for n in include])}" that have diverged or are absent in "{", ".join([eggnog.get_species_name_by_id(n, df_species) for n in exclude])}"'
    )

# create temporary .csv of homologs for further use
homologs_df.reset_index().to_csv("temp/1_A_homologs_.tsv", sep="\t", index=False)

###----------------------------------------------------------------------
### 1) B) extract unique protein IDs from homologs file and

print("Extracting unique protein IDs from homologs ...")

# extract unique protein IDs from homologs dataframe
unique_protein_ids = (
    homologs_df["protein_id"].str.split(",").explode().str.strip().unique()
)

print("Writing unique protein IDs to 1_B_unique_protein_IDs.txt file in /results ...")

with open("results/1_B_unique_protein_IDs.txt", "w") as f:
    f.write("\n".join(unique_protein_ids))

###---------------------------------------
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
    functional_categories_in_homologs.apply(list)
    .explode()
    .value_counts()
    .reset_index(name="count")
    .rename(columns={"index": "category_code"})
)

category_counts_df = category_counts_df.merge(
    df_functional_categories_description, on="category_code", how="left"
)

print("Exporting result table to results/1_C_functional_category_counts.csv ...")

category_counts_df.to_csv("results/1_C_functional_category_counts.csv", index=False)
