## applying partition model
iqtree2 -s /fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree/sanitized_sequences.fasta -p /fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree/sanitized_sequences.fasta.splits.nex -B 1000 -T AUTO

## Choosing the best partitioning scheme
iqtree2 -s turtle.fa -p turtle.nex -B 1000 -T AUTO -m MFP+MERGE -rcluster 10 --prefix turtle.merge

## Tree topology tests
### First, concatenate the trees constructed by single and partition models into one file:
cat turtle.fa.treefile turtle.nex.treefile >turtle.trees
### pass this file into IQ-TREE via -z option:
iqtree2 -s turtle.fa -p turtle.merge.best_scheme.nex -z turtle.trees -zb 10000 -au -n 0 --prefix turtle.test