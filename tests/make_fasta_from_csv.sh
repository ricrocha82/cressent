awk -F , 'NR > 1 {print ">" $1 "|" $3 "\n" $2}' /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/ORF_seq.csv > /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/reps.fa


# searching fro walker a motif
seqkit locate --id-ncbi -i -r -p "[GA].{4}GK[TS]" /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps_align_trimm.fa > /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_location.txt

# split the rep seqs






# seqkit locate -i -r -p "[GA].{4}GK[TS]" /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps.fa \
#     --gtf > /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps_walkera.gtf

# seqkit subseq --only-flank --gtf  \
#                 /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps_walkera.gtf \
#                 /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps.fa \
#                 > /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps_split.fa
