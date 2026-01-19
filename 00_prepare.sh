#!/bin/bash

mkdir -p results
mkdir -p doc
mkdir -p temp
mkdir -p data
mkdir -p data/zips

wget -P data/zips/ http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_members.tsv.gz 
wget -P data/zips/ http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_annotations.tsv.gz
wget -P data/ http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv
wget -P data/ http://eggnog5.embl.de/download/eggnog_4.5/eggnog4.functional_categories.txt

for f in data/zips/*.gz; do
    base_name=$(basename "$f") # Get the filename of zips without the path    
    gunzip -c "$f" > data/"${base_name%.gz}" # unzip and save to data folder
    echo "Zip-file $base_name was extracted to data folder." # echo info message to terminal
done