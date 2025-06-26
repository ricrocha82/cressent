#!/bin/bash

# I used this script because I forgot to use the --keep_names

BASE_DIR="/fs/project/PAS1117/ricardo/cressent/DB/caps_cluster"
cd "$BASE_DIR" || exit

for dir in */ ; do
    family=$(basename "$dir")
    tsv="${dir}${family}_name_table.tsv"
    fasta="${dir}renamed_${family}.fa"
    output="${dir}${family}.fa"

    if [[ -f "$tsv" && -f "$fasta" ]]; then
        echo "Processing $family..."

        python3 - <<EOF
import csv

# Build mapping
mapping = {}
with open("$tsv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sanitized = row["Sanitized_Name"].strip()
        original = row["Original_Name"].strip()
        mapping[sanitized] = original

# Rename headers in FASTA
with open("$fasta", "r") as fin, open("$output", "w") as fout:
    for line in fin:
        if line.startswith(">"):
            header = line[1:].strip()
            new_header = mapping.get(header, header)
            fout.write(f">{new_header}\n")
        else:
            fout.write(line)
EOF

    else
        echo "Skipping $family: missing .tsv or .fa file"
    fi
done
