#!/bin/bash

#SBATCH --account=PAS1117
#SBATCH --job-name=cluster_caps.sh
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --error=/fs/project/PAS1117/ricardo/cressent/DB/scripts/jobs/%x_%j.err
#SBATCH --output=/fs/project/PAS1117/ricardo/cressent/DB/scripts/jobs/%x_%j.out

source ~/miniconda3/bin/activate
conda activate test-cressent

INPUT_DIR="/fs/project/PAS1117/ricardo/cressent/DB/caps"
OUTPUT_DIR="/fs/project/PAS1117/ricardo/cressent/DB/caps_cluster"
SCRIPT="/fs/project/PAS1117/ricardo/cressent/cressent/cli.py"
THREADS=40

for file in "$INPUT_DIR"/*.fa; do
    filename=$(basename "$file")                # e.g., Genomoviridae.fa
    base="${filename%.fa}"                      # e.g., Genomoviridae
    input_path="$file"
    output_path="$OUTPUT_DIR/$base"

    echo "Processing $filename..."
    python "$SCRIPT" cluster -t "$THREADS" -i "$input_path" -o "$output_path" --keep_names
done

echo "job is done"