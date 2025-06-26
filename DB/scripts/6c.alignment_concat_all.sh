
# source ~/miniconda3/bin/activate
# conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/familycopy.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db

output_dir="$dir_to_download/aligned_db/all/"

mkdir -p "$output_dir"



# REPS
input_dir_reps="$dir_to_download/concat/reps/"
# Loop through each project name from the family file
for fasta_file in "$input_dir_reps"/*.fa; do
    # Define the output file path
    cat $fasta_file >> "$output_dir/all_rep.fa"

    # align and trim
    output_file="$output_dir/aligned_all_rep.fa"  
    output_file_trim="$output_dir/aligned_all_rep_trim.fa"

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"
    trimal -keepheader -in "$output_file" -out "$output_file_trim" -gt 0.15

done

# CAPS
input_dir_cap="$dir_to_download/concat/caps/"
# Loop through each project name from the family file
for fasta_file in "$input_dir_cap"/*.fa; do
    # Define the output file path
    cat $fasta_file >> "$output_dir/all_cap.fa"

    # align and trim
    output_file="$output_dir/aligned_all_cap.fa"  
    output_file_trim="$output_dir/aligned_all_cap_trim.fa"

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"
    trimal -keepheader -in "$output_file" -out "$output_file_trim" -gt 0.15

done
