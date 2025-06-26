
# source ~/miniconda3/bin/activate
# conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA/data/family.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db
# input_dir="$dir_to_download/annotated/caps/"
# output_dir="$dir_to_download/aligned/caps/"
input_dir="$dir_to_download/annotated/reps/"
output_dir="$dir_to_download/aligned/reps/"

mkdir -p "$output_dir"

# Loop through each project name from the family file
for fasta_file in "$input_dir"/*.fa; do
    # Extract the base name of the file (without directory and extension)
    base_name=$(basename "$fasta_file" .fa)

    # Define the output file path
    output_file="$output_dir/${base_name}_aligned.fa"

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file"

done
