source ~/miniconda3/bin/activate
conda activate biopython

family_file=/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/names_to_query_ssdna.csv
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/families_aa/

for project in `awk '{print $1}' "$family_file"`
do
    echo 'download files from' ${project}
    mkdir -p $dir_to_download/$project
    # search_term="($project [Organism] OR $project [All Fields]) AND capsid protein [Title]"
    search_term="($project [Organism] OR $project [All Fields])"
    echo $search_term
    python /fs/project/PAS1117/ricardo/ssDNA_test/scripts/download_from_ncbi_family.py -e 'pavan.4@osu.edu' -t "$search_term" >> $dir_to_download/${project}.fa
    
    # python /fs/project/PAS1117/ricardo/ssDNA/scripts/download_from_ncbi_family_na.py -e 'pavan.4@osu.edu' -t "$search_term" >> $dir_to_download/$project/${project}_na.fa
    # sed -i '/^$/d' $dir_to_download/$project/${project}.fa
done

# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa
