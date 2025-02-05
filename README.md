# ssDNA tool
a modular tool to help researchers to automatically annotate ssDNA contigs

# Install
Use requirements.txt

# Test data
the test data is a fasta file from [Kazlauskas et al 2018](https://www.mdpi.com/1999-4915/10/4/187)

# Run the codes
- Steps:

  1. align.py (MAFFT -> Trimal) (./ssDNA_tool/ssDNA_annotator/modules/align.py)
  2. split sequence based on determined pattern/motif (seqkit and custom script) (./ssDNA_tool/ssDNA_annotator/modules/pattern_find_split.py)
  3. make phylogenetic tree (IQTREE) (./ssDNA_tool/ssDNA_annotator/modules/build_tree.py)
  4. tanglegram (R) (./ssDNA_tool/ssDNA_annotator/modules/tanglegram.py)
  5. Sequence logo (R) (./ssDNA_tool/ssDNA_annotator/modules/seq_logo.py)
  6. Cluster sequences (./ssDNA_tool/ssDNA_annotator/modules/cluster.py - updated to Python3)

Code examples to run each file (arguments, etc) is at the end of each python script

# DB with all protein sequences
[Families](https://github.com/ricrocha82/ssDNA_tool/blob/main/ssDNA_annotator/test_data/family.txt) were collected from the ICTV website

Only the Pleolipoviridae Family (contains both dsDNA and ssDNA) were filtered by Genus names:
- Alphapleolipovirus finnoniense
- Alphapleolipovirus huluense
- Alphapleolipovirus samutsakhonense
- Alphapleolipovirus thailandense

Sequences that were not found in Refseq were directly download (if available) using the Acession number extracted from the reports found on [ictv report website](https://ictv.global/report):
- Adamaviridae
- Anicreviridae
- Draupnirviridae
- Endolinaviridae (not found)
- Gandrviridae
- Geplanaviridae (not found)
- Kanorauviridae (not found)
- Kirkoviridae
- Mahapunaviridae (not found)
- Ouroboviridae
- Pecoviridae

All sequences are here: "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/families_aa"
All the families are here: "/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/family.txt"

## processing DB (for global comparison)
- clustered at 95% amino acid identity > 90% (`CD-HIT`)
- clustered using `MCL` algorithm - inflation rate at 1.5
- Clusters with proteins >= 10 annotated Rep and Capsid will be aligned using `MAFFT` with the auto parameter.


