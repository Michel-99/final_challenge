import eggnog_library as eggnog
import os

# Dataframe setup
df_annotations = eggnog.dataframe_setup_annotations()
df_members = eggnog.dataframe_setup_members()
df_species = eggnog.dataframe_setup_taxid_info()
df_functional_categories = eggnog.dataframe_setup_functional_categories()

###------
### 1) A) homologuous genes in humans and chimp but not mouse

# get species ids
human_id = eggnog.get_species_id_by_name("Homo sapiens", df_species)
chimp_id = eggnog.get_species_id_by_name("Pan troglodytes", df_species)
mouse_id = eggnog.get_species_id_by_name("Mus musculus", df_species)

# set which species should be included and excluded
include = [human_id, chimp_id]
exclude = [mouse_id]

# look for homologs in included species (exclude excluded species)
homologs = eggnog.filter_by_ids(
    df_members, "species_taxid_containing_protein", include, exclude
)

# output to .txt file
with open("results/1_A.txt", "w") as f:
    f.write(
        f"The species {include} share {len(homologs)} homolog genes that {exclude} don´t have"
    )

homologs.to_csv("results/1_A_homologs_.tsv", sep="\t", index=False)
