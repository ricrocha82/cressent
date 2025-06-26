
source ~/miniconda3/bin/activate
conda activate cd-hit

# make a new family.txt 
# for file in *.fa; do basename "$file" .fa; done > filenames.txt

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/family.txt
dir_db=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/
dir_aa="$dir_db/families_aa"

for project in `awk '{print $1}' "$family_file"`
do
    echo 'Cluster sequences from' ${project}
    mkdir -p $dir_db/cd_hit/
    in=$dir_aa/${project}.fa
    out=$dir_db/cd_hit/${project}.fa

    # cluster at 95% amino acid identity (-c) >90% of the shortest sequence (-aS)
    cd-hit -i $in -o $out -c 0.95 -aS 0.9 -M 16000 -T 8

done

# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa
