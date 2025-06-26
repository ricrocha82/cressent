#!/bin/bash
# source ~/miniconda3/bin/activate
# conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/familycopy.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db

output_dir="$dir_to_download/aligned_db/all/"

mkdir -p "$output_dir"

# REPS
# Loop through each project name from the family file
for project in `awk '{print $1}' "$family_file"`
do
    # Define the output file path
    fasta_file="$input_dir/${project}_concat_files.fa"
    output_file="$output_dir/${project}_aligned_rep.fa"  
    # output_file="$output_dir/${project}_aligned_cap.fa"  

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"

done
trimal -keepheader -in {aligned_fasta} -out {trimmed_file} -gt 0.15

# REPS
input_dir_reps="$dir_to_download/concat/reps/"
# Loop through each project name from the family file
for fasta_file in "$input_dir_reps"/*.fa; do
    # Define the output file path
    cat $fasta_file >> "$output_dir/all_rep.fa"
    output_file="$output_dir/aligned_all_rep.fa"  

    # cat $fasta_file >> "$output_dir/${project}_all_cap.fa"
    # output_file="$output_dir/${project}_aligned_all_cap.fa"  

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"

done

# CAPS
input_dir_cap="$dir_to_download/concat/caps/"
# Loop through each project name from the family file
for fasta_file in "$input_dir_cap"/*.fa; do
    # Define the output file path
    cat $fasta_file >> "$output_dir/all_cap.fa"
    output_file="$output_dir/aligned_all_cap.fa"  

    # cat $fasta_file >> "$output_dir/${project}_all_cap.fa"
    # output_file="$output_dir/${project}_aligned_all_cap.fa"  

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"

done
