# Metazoan Gene Conservation Analysis
## Function
This project analyzes gene conservation patterns across metazoan species using EggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) database. 
The analysis focuses on:
- Identifying primate-specific genes
- Analyzing gene loss in rodent lineages
- Discovering universally conserved animal genes

## Tools
Programming language: Python 3.9.12
Third party libraries: pandas, os, sys, csv
Own created libraries: eggnog_library

## How to run
1. execute ```bash runall.sh```
or 
1. execute ```bash runall.sh --skip-analysis```
2. execute ```python main.py```
## Structure
runall.sh
1. create directory structure 
   * /data
   * /results
  
2. download eggnog data for metazoans (ID 33208) into /data directory

main.py
3. set up result .txt file structure in /results directory and assign variables to result files for better usability
4. create pandas dataframes from previously downloaded third party eggnog data 
5. look for homologs genes in humans and chimps but not mice
   * Get species id from name (string)
   * define which species IDs should be included and excluded
   * filter database for defined include/exclude
   * output to result-file
6. extract protein IDs for found homologs from previous step
   * extract from pandas dataframe and convert to single line output string (using .strip(), .unique(), etc.)
7. analyze functional categories of homologous genes 
   * merge homologs with annotations on orthologous_group_id
   * explode functional categories (genes can have multiple)
   * count occurrences and merge with category descriptions
   * export to CSV with readable category names
8. identify OGs unique to Human and Chimp 
   * filter for OGs where num_of_species == 2, using the already for human and chimp filtered dataset
   * isolate genes present only in this species pair
9. identify primate-specific OGs (1E)
   * use filter_by_species_names() with PRIMATES list
   * find OGs conserved across all primates but absent in other lineages


## Troubleshooting
#Import errors
- Ensure all dependencies are installed: `pip install pandas`
- Check Python version: `python --version` (should be 3.9+)
#Empty result files
- Verify data downloaded correctly in `/data` directory
- Check that 00_prepare.sh completed without errors

    ## Data Source
This analysis uses the EggNOG 5.0 database:
- Database: [EggNOG](http://eggnog5.embl.de/)
- Taxonomic level: Metazoa (TaxID: 33208)
- Citation: Huerta-Cepas et al. (2019) eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses.
Nucleic Acids Res. 47:D309-D314

## Results
QUESTION 1: PRIMATE-SPECIFIC GENE ANALYSIS
1A) Homologs in Human & Chimp but NOT in Mouse: 1073 OGs
1B) Unique protein IDs: 52188
1C) Functional categories found: 23
    Top 3 categories:
      - S: 534 genes
      - K: 162 genes
      - T: 131 genes
1D) OGs found ONLY in Human and Chimp: 21 OGs
1E) Primate-specific OGs: 167 OGs

QUESTION 2: LINEAGE ANALYSIS
Core vertebrate OGs (Primates + Chicken + Fish): 9861
Lost in BOTH Mouse and Rat: 108
Lost ONLY in Mouse: 63
Lost ONLY in Rat: 138

QUESTION 3: UNIVERSAL ANIMAL GENES
Total unique species in dataset: 161
Universal OGs (99%+): 694

# group members
Josef Birnöcker: 1)a-c)
Clara Pernold: 2 & 3
Michel Zwicker: 1)b-e)
All of the individuals mentioned contributed equally to the project, and the codes and ideas were regularly discussed in team meetings.

## References
1. EggNOG 5.0 Database - http://eggnog5.embl.de/
2. https://biopython.org/docs/latest/Tutorial/index.html
3. https://stackoverflow.com/questions/11707586/how-do-i-expand-the-output-display-to-see-more-columns-of-a-pandas-dataframe
4. https://stackoverflow.com/questions/642154/how-do-i-convert-all-strings-in-a-list-of-lists-to-integers
5. https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
6. https://www.geeksforgeeks.org/python/different-types-of-joins-in-pandas/
7. https://stackoverflow.com/questions/19960077/how-to-filter-pandas-dataframe-using-in-and-not-in-like-in-sql
8. https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html
9. https://stackoverflow.com/questions/34962104/how-can-i-use-the-apply-function-for-a-single-column
10. https://www.geeksforgeeks.org/python/python-set-operations-union-intersection-difference-symmetric-difference/
11. https://docs.python.org/3/tutorial/datastructures.html#sets
12. https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.explode.html
13. https://pandas.pydata.org/docs/reference/api/pandas.Series.str.split.html
14. https://stackoverflow.com/questions/26147180/convert-row-to-column-header-for-pandas-dataframe
15. https://stackoverflow.com/questions/20490274/how-to-reset-index-in-a-pandas-dataframe
16. https://www.w3schools.com/python/pandas/default.asp
17. https://realpython.com/pandas-dataframe/
18. https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe
