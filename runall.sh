#!/bin/bash

set -e  

# Parse command line arguments
RUN_ANALYSIS=true
if [[ "$1" == "--skip-analysis" ]]; then
    RUN_ANALYSIS=false
fi

echo "Creating directory structure..."
mkdir -p results
mkdir -p doc
mkdir -p data
mkdir -p data/zips

echo "Downloading EggNOG data files..."

wget -P data/zips/ http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_members.tsv.gz 
wget -P data/zips/ http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_annotations.tsv.gz
wget -P data/ http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv
wget -P data/ http://eggnog5.embl.de/download/eggnog_4.5/eggnog4.functional_categories.txt

echo "Extracting compressed files..."
for f in data/zips/*.gz; do
    base_name=$(basename "$f")
    gunzip -c "$f" > data/"${base_name%.gz}"
    echo "Extracted: $base_name"
done

echo ""
echo "Setup complete!"

if [ "$RUN_ANALYSIS" = true ]; then
    echo ""
    echo "Running analysis..."
    python 1-3_runall.py
    echo ""
    echo "All done! Check results/ folder for output files."
else
    echo ""
    echo "Skipping analysis. Run manually with: python main.py"
fi

