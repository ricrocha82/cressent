source ~/miniconda3/bin/activate
conda activate sra_tools

# this code has a weird behaviour naming directories with an ? at the end (changed to a python script)

family_file=/fs/project/PAS1117/ricardo/ssDNA_test/test_data/accession_not_found_refseq.csv
dir_to_download=/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/

awk -F',' 'NR>1 {print $1, $2}' "$family_file" > "$dir_to_download/temp_file.txt"

while read -r project family
do
    echo "Processing: Accession = ${project}, Family = ${family}"

    search_term="${project}"
    mkdir -p $dir_to_download/$family
    # outputs
    out_1=$dir_to_download/${family}/${project}.txt
    out_2=$dir_to_download/${family}/${family}.fa

    # Redirect stdin to /dev/null to avoid interference
    ( esearch -db nuccore -query "$search_term" | efetch -format ft | awk '/protein_id/ {print $2}' > "$out_1" ) < /dev/null

    for protein_id in $(awk '{print $1}' $out_1); do
    echo $protein_id
    (esearch -db protein -query "$protein_id" | efetch -format fasta >> $out_2) < /dev/null
    done

    rm $out_1

done < "$dir_to_download/temp_file.txt"

rm "$dir_to_download/temp_file.txt"


# get a coding region (CDS) from a complete genome
# esearch -db nuccore -query 'KM821764.1' | efetch -format fasta_cds_na > /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa

# python /fs/project/PAS1117/ricardo/ssDNA/scripts/test1.py -seq /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa -out /fs/project/PAS1117/ricardo/ssDNA/data/CDS_transl.fa

# python /fs/project/PAS1117/ricardo/ssDNA/scripts/translate_seq.py -i /fs/project/PAS1117/ricardo/ssDNA/data/CDS.fa -o /fs/project/PAS1117/ricardo/ssDNA/data/CDS_transl.fa

# esearch -db nuccore -query KM821764.1 | efetch -format ft | awk '/protein_id/ {print $2}' > /fs/project/PAS1117/ricardo/ssDNA/data/protein_ids.txt

# for protein_id in $(awk '{print $1}' /fs/project/PAS1117/ricardo/ssDNA/data/protein_ids.txt); do
#     echo $protein_id
#     esearch -db protein -query "$protein_id" | efetch -format fasta >> /fs/project/PAS1117/ricardo/ssDNA/data/protein_ids.fasta
# done

# cat gb*.fasta > all_proteins.fasta

# esearch -db nuccore -query KM821764.1 | elink -target protein | esummary


# esearch -db nuccore -query KM821764.1 | elink -target protein |  esummary | xtract -pattern DocumentSummary -if Title -contains "CDS" -element Extr



# esearch -db nuccore -query 'Anelloviridae [Organism] OR Anelloviridae [All Fields]' | efetch -format fasta >> /fs/scratch/Sullivan_Lab/Ricardo/Anelloviridae/Anelloviridae_na.fa