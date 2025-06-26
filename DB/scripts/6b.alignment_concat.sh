
# source ~/miniconda3/bin/activate
# conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/familycopy.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db


# REP
input_dir="$dir_to_download/concat/reps/"
output_dir="$dir_to_download/aligned_db/reps/"

mkdir -p "$output_dir"

# Loop through each project name from the family file
for project in `awk '{print $1}' "$family_file"`
do
    # Define the output file path
    fasta_file="$input_dir/${project}_concat_files.fa"
    output_file="$output_dir/${project}_aligned_rep.fa"  
    output_file_trim="$output_dir/${project}_aligned_rep_trim.fa" 

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"
    trimal -keepheader -in "$output_file" -out "$output_file_trim" -gt 0.15

done



# CAP
input_dir_cap="$dir_to_download/concat/caps/"
output_dir_cap="$dir_to_download/aligned_db/caps/"
mkdir -p "$output_dir_cap"

# Loop through each project name from the family file
for project in `awk '{print $1}' "$family_file"`
do
    # Define the output file path
    fasta_file="$input_dir_cap/${project}_concat_files.fa" 
    output_file="$output_dir_cap/${project}_aligned_cap.fa"  
    output_file_trim="$output_dir_cap/${project}_aligned_cap_trim.fa" 

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"
    trimal -keepheader -in "$output_file" -out "$output_file_trim" -gt 0.15

done
