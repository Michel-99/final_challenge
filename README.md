# Metazoan Gene Conservation Analysis
## Function
Compare conserved genes in metazoans between primates and other mammals


## Tools
Programming language: Python 3.9.12
Third party libraries: pandas, os, sys, csv
Own created libraries: eggnog_library

## How to run
1. execute bash 00_prepare.sh
2. execute runall.py

## Structure
00_prepare.ssh
1. create directory structure 
   * /data
   * /results
   * temp
2. download eggnog data for metazoans (ID 33208) into /data directory

runall.py
3. set up result .txt file structure in /results directory and assign variables to result files for better usability
4. create pandas dataframes from previously downloaded third party eggnog data 
5. look for homologs genes in humans and chimps but not mice
   * Get species id from name (string)
   * define which species IDs should be included and excluded
   * filter database for defined include/exclude
   * output to result-file
6. extract protein IDs for found homologs from previous step
   * extract from pandas dataframe and convert to single line output string (using .strip(), .unique(), etc.)




## Results


# group members
Josef Birnöcker: 1)a)
Clara Pernold: 2 & 3
Michel Zwicker: 1)b-e)
