
source ~/miniconda3/bin/activate
conda activate diamond

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/family.txt
dir_db=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/


for project in `awk '{print $1}' "$family_file"`
do
    echo 'Cluster sequences from' ${project}
    mkdir -p $dir_db/diamond/
    in=$dir_db/cd_hit/${project}.fa
    out_1=$dir_db/diamond/${project}
    out_2=$dir_db/diamond/${project}.txt

    # compare sequences all-vs-all evalue threshol=0.00001
    diamond makedb --in $in -d $out_1 --threads 24

    diamond blastp --db $out_1 \
                    --query $in \
                    --out $out_2 \
                    --outfmt 6 --threads 24 --evalue 0.00001 \
                    --masking 0 --sensitive 

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