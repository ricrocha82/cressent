
source ~/miniconda3/bin/activate
conda activate mcl

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/names_to_query_ssdna.csv
dir_db=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db

for project in `awk '{print $1}' "$family_file"`
do
    echo 'Cluster sequences from' ${project}
    mkdir -p $dir_db/mcl/
    in=$dir_db/diamond/${project}.txt
    out_1=$dir_db/diamond/${project}_out
    out_2=$dir_db/mcl/${project}_mcl_evalue

    # preapre data for MCL input
    # e-value as our similarity measure
    cut -f1,2,11 $in > $out_1

    # run MCL
    # The output is then a file where each line is a cluster of tab-separated labels
    mcl $out_1 --abc -I 1.5 -te 24 -o $out_2 --abc-neg-log10

done

# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa

# diamond makedb --in /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/cd_hit/Bidnaviridae.fa \
#                     -d /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/diamond/Bidnaviridae --threads 8


# diamond blastp --masking 0 --sensitive --db /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/diamond/Bidnaviridae --out /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/diamond/Bidnaviridae.txt --outfmt 6 --threads 8 --evalue 0.00001


# diamond blastp --db /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/diamond/Bidnaviridae \
#                     --query /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/cd_hit/Bidnaviridae.fa \
#                     --out /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/diamond/Bidnaviridae.txt \
#                     --outfmt 6 --threads 24 --evalue 0.00001 \
#                     --masking 0 --sensitive 

