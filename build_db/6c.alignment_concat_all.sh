#! usr/bin/env bash


source ~/miniconda3/bin/activate
conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/names_to_query_ssdna.csv
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db

output_dir="$dir_to_download/aligned_db/all/"

mkdir -p "$output_dir"

# Concatenate REPS
input_dir_reps="$dir_to_download/concat/reps/"
input_dir_caps="$dir_to_download/concat/caps/"
# Loop through each project name from the family file

# Concatenate all Reps into a FASTA file
# cat "$input_dir_reps"/*.fa >> "$output_dir/all_rep.fa"
# # Concatenate all Caps into a FASTA file
# cat "$input_dir_caps"/*.fa >> "$output_dir/all_cap.fa"

## RUN MAFFT and TRIMAL in all the Rep and Cap sequences
# align and trim
fasta_file="$output_dir/all_rep.fa"
output_file="$output_dir/aligned_all_rep.fa"  
output_file_trim="$output_dir/aligned_all_rep_trim.fa"
# # Run MAFFT alignment
mafft --auto "$fasta_file" > "$output_file"
trimal -keepheader -in "$output_file" -out "$output_file_trim" 


fasta_file="$output_dir/all_cap.fa"  
output_file="$output_dir/aligned_all_cap.fa"  
output_file_trim="$output_dir/aligned_all_cap_trim.fa"
# # Run MAFFT alignment
mafft --auto "$fasta_file" > "$output_file"
trimal -keepheader -in "$output_file" -out "$output_file_trim" -gt 0.15