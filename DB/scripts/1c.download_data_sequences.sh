source ~/miniconda3/bin/activate
conda activate sra_tools

family_file=/fs/project/PAS1117/ricardo/ssDNA/data/family.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/families_nc/

for project in `awk '{print $1}' "$family_file"`
do
    echo 'download files from' ${project}
    mkdir -p $dir_to_download/$project
    # search_term="($project [Organism] OR $project [All Fields]) AND capsid protein [Title]"
    search_term="($project [Organism] OR $project [All Fields])"
    echo $search_term
    esearch -db nuccore -query $search_term | efetch -format fasta >> $dir_to_download/$project/${project}_nc.fasta
done

# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa
