#!/bin/bash

SOURCE_BASE="/fs/project/PAS1117/ricardo/cressent/DB/caps_cluster"
DEST_DIR="/fs/project/PAS1117/ricardo/cressent/DB/caps_new"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each subdirectory
for dir in "$SOURCE_BASE"/*/; do
    family=$(basename "$dir")
    file="${dir}${family}.fa"
    
    if [[ -f "$file" ]]; then
        echo "Moving $file → $DEST_DIR/$family.fa"
        mv "$file" "$DEST_DIR/$family.fa"
    else
        echo "Skipping $family: $family.fa not found"
    fi
done