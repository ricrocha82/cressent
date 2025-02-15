
#!/usr/bin/env bash

source ~/miniconda3/bin/activate
conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/names_to_query_ssdna.csv
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db


## CAPS
input_dir_caps="$dir_to_download/annotated/caps/"
output_dir_caps="$dir_to_download/concat/caps/"
output_dir_aligned_caps="$dir_to_download/aligned_db/caps/"

mkdir -p "$output_dir_caps"
mkdir -p "$output_dir_aligned_caps"

for project in `awk '{print $1}' "$family_file"`
do

    echo -e "\n Cap Family $project \n"

    mkdir -p "$dir_to_download/temp/${project}/"

    # Concatenate all project-specific FASTA files
    cat "$input_dir_caps"/${project}*.fa > "$dir_to_download/temp/${project}/${project}.fa"
    # add family name at the end of the sequence name
    python /fs/project/PAS1117/ricardo/ssDNA_tool/build_db/6b.add_family_name.py -i "$dir_to_download/temp/${project}/" -o "$output_dir_caps"

    rm -r "$dir_to_download/temp/${project}/"

    # Define the output file path
    fasta_file="$output_dir_caps/${project}.fa"
    output_file_caps="$output_dir_aligned_caps/${project}_aligned_cap.fa"  
    output_file_trim_caps="$output_dir_aligned_caps/${project}_aligned_cap_trim.fa" 

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file_caps"
    trimal -keepheader -in "$output_file_caps" -out "$output_file_trim_caps" -gt 0.15
done

## REPS
input_dir_reps="$dir_to_download/annotated/reps/"
output_dir_reps="$dir_to_download/concat/reps/"
output_dir_aligned_reps="$dir_to_download/aligned_db/reps/"

mkdir -p "$output_dir_reps"
mkdir -p "$output_dir_aligned_reps"

for project in `awk '{print $1}' "$family_file"`
do
    echo -e "\n Rep Family $project \n"

    mkdir -p "$dir_to_download/temp/${project}/"

    # Concatenate all project-specific FASTA files
    cat "$input_dir_reps"/${project}*.fa > "$dir_to_download/temp/${project}/${project}.fa"

    # add family name at the end of the sequence name
    python /fs/project/PAS1117/ricardo/ssDNA_tool/build_db/6b.add_family_name.py -i "$dir_to_download/temp/${project}/" -o "$output_dir_reps"

    rm -r "$dir_to_download/temp/${project}/"

    # Define the input and output file paths
    fasta_file="$output_dir_reps/${project}.fa"
    output_file_reps="$output_dir_aligned_reps/${project}_aligned_rep.fa"  
    output_file_trim_reps="$output_dir_aligned_reps/${project}_aligned_rep_trim.fa" 

    # # Run MAFFT alignment
    mafft --auto "$fasta_file" > "$output_file_reps"
    trimal -keepheader -in "$output_file_reps" -out "$output_file_trim_reps" -gt 0.15
done