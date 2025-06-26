source ~/miniconda3/bin/activate
conda activate sra_tools

family_file=/fs/project/PAS1117/ricardo/ssDNA/data/family.txt
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/complete_genome/

for project in `awk '{print $1}' "$family_file"`
do
    echo 'download files from' ${project}
    mkdir -p $dir_to_download/$project
    # search_term="($project [Organism] OR $project [All Fields]) AND capsid protein [Title]"
    search_term="($project [Organism] OR $project [All Fields])"
    echo $search_term
    esearch -db assembly -query $search_term | esummary \
    | xtract -pattern DocumentSummary \
        -element AssemblyAccession \
        -element AssemblyName \
        -element SpeciesName \
        -element FtpPath_RefSeq > $dir_to_download/$project/${project}.txt
    sed -i '/^GCF_/!d' $dir_to_download/$project/${project}.txt
    for protein_id in $(awk '{print $1}' $dir_to_download/$project/${project}.txt); do
        echo $protein_id
        esearch -db assembly -query $protein_id | elink -target nuccore | efetch -format fasta >> $dir_to_download/$project/${project}_complete_genome.fasta
    done
    # python /fs/project/PAS1117/ricardo/ssDNA/scripts/download_from_ncbi_family_na.py -e 'pavan.4@osu.edu' -t "$search_term" >> $dir_to_download/$project/${project}_na.fa
    # sed -i '/^$/d' $dir_to_download/$project/${project}.fa
done

# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa

# esearch -db assembly -query 'Anelloviridae' | esummary \
#     | xtract -pattern DocumentSummary \
#         -element AssemblyAccession \
#         -element AssemblyName \
#         -element SpeciesName \
#         -element FtpPath_RefSeq > /fs/project/PAS1117/ricardo/ssDNA/data/test.txt

# sed -i '/^GCF_/!d' /fs/project/PAS1117/ricardo/ssDNA/data/test.txt

#     for protein_id in $(awk '{print $1}' /fs/project/PAS1117/ricardo/ssDNA/data/test.txt); do
#     echo $protein_id
#     esearch -db assembly -query $protein_id | elink -target nuccore | efetch -format fasta >> /fs/project/PAS1117/ricardo/ssDNA/data/${project}_complete_genome.fasta
#     done


