# source ~/miniconda3/bin/activate
# conda activate biopython

from Bio import SeqIO
import os
import glob

# Directory containing the FASTA files
input_dir = "path/to/fasta/files"  # Change this to your directory
output_dir = "path/to/output/files"  # Change this to your output directory

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process all FASTA files in the directory
for fasta_file in glob.glob(os.path.join(input_dir, "*.fasta")):
    family_name = os.path.basename(fasta_file).replace(".fasta", "")  # Extract family name from filename
    output_file = os.path.join(output_dir, os.path.basename(fasta_file))  # Keep same filename in output

    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.id = f"{record.id}|{family_name}"  # Append family name to header
            record.description = ""  # Remove description to avoid duplication
            SeqIO.write(record, out_f, "fasta")

print("Processing complete! FASTA files with updated headers are saved in:", output_dir)