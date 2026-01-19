#!/bin/bash

mkdir -p results
mkdir -p doc
mkdir -p temp
mkdir -p data
cd data

wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_members.tsv.gz
wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_annotations.tsv.gz
wget http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv
wget http://eggnog5.embl.de/download/eggnog_4.5/eggnog4.functional_categories.txt