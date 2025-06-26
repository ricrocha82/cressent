# Methods to generate the database of the ssDNA tool

1 - Protein sequences will be retrieved at Family levels using the `entrez API` with the following query: {family_name} [Organism] OR {family_name} [All Fields]
(/fs/project/PAS1117/ricardo/ssDNA/scripts/1.download_data.sh)

1.a - CDS will be retrived from complete genomes at Genus kevel (links from ICTV: https://ictv.global/report) then translated into protein sequences using `SeqIO`
(/fs/project/PAS1117/ricardo/ssDNA/scripts/1.download_data_CDS.sh) 

2 - To reduce redundancy, sequences will be clustered at 95% amino acid identity > 90% (`CD-HIT`)
(/fs/project/PAS1117/ricardo/ssDNA/scripts/2.cluster_cdhit.sh)

3 - Resultant sequences will be compared against DIAMOND (e-value 1e-5) using `BALSTp` all-vs-all and against PDB, UniRef90 and Swiss-Prot
(/fs/project/PAS1117/ricardo/ssDNA/scripts/3.blastp.sh)

4 - The result will be clustered using `MCL` algorithm - inflation rate at 1.5
4.1 - Sequences were separated in three groups: 
    - Capsid and VP-related proteins
    - Rep proteins
    - unannotated ORFs
(/fs/project/PAS1117/ricardo/ssDNA/scripts/4.cluster_MCL.sh)
(/fs/project/PAS1117/ricardo/ssDNA/scripts/5a.cluster_protein_description.R)
(/fs/project/PAS1117/ricardo/ssDNA/scripts/5b.cluster_protein_curation.R)

5 - use cluster module to reduce redundancy in each family