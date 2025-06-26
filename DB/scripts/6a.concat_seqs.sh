
# source ~/miniconda3/bin/activate
# conda activate mafft

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/familycopy.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db
input_dir="$dir_to_download/annotated/caps/"
output_dir="$dir_to_download/concat/caps/"
# input_dir="$dir_to_download/annotated/reps/"
# output_dir="$dir_to_download/concat/reps/"

mkdir -p "$output_dir"

for project in `awk '{print $1}' "$family_file"`
do
    for fasta_file in "$input_dir"/${project}*.fa; do

    # # Run MAFFT alignment
    cat $fasta_file >> "$output_dir/${project}_concat_files.fa"

done

done