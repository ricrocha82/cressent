# ssDNA tool
a modular tool to help researchers to automatically annotate ssDNA contigs

# Install
Use requirements.txt

# Test data
the test data is a fasta file from [Kazlauskas et al 2018](https://www.mdpi.com/1999-4915/10/4/187)

# Run the codes
- Steps:

    1- align.py (MAFFT -> Trimal) (/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py)
    2- split sequence based on determined pattern/motif (seqkit and custom script) (/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/pattern_find_split.py)
    3 - make phylogenetic tree (IQTREE) (/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py)
    4 - tanglegram (R) (/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/tanglegram.py)
    5 - Sequence logo (R) (/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/seq_logo.py)

The code to run each file (arguments, etc) is at the end of each python script

