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
  3. make phylogenetic tree (using [`IQTREE`](http://www.iqtree.org/)) (./ssDNA_tool/ssDNA_annotator/modules/build_tree.py)
  4. tanglegram (R) (./ssDNA_tool/ssDNA_annotator/modules/tanglegram.py)
  5. Sequence logo (R) (./ssDNA_tool/ssDNA_annotator/modules/seq_logo.py)
  6. Cluster sequences (./ssDNA_tool/ssDNA_annotator/modules/cluster.py)

Code examples to run each file (arguments, etc) is at the end of each python script

# References
## for phylogenetic tree analysis using R
- https://yulab-smu.top/treedata-book




---
# ssDNA Annotator Pipeline Example

This pipeline demonstrates how to process metagenomic sequences and build phylogenetic trees using the ssDNA_annotator tool. The workflow includes sequence clustering (dereplication), alignment, tree building, and visualization.


## Input Files

- **Metagenomic Sequences FASTA**  
/path/to/my_sequence.fa

- **Database for Phylogenetic Analysis**  
path/to/database

- here [build_db](https://github.com/ricrocha82/ssDNA_tool/tree/main/build_db) is the pipeline used to build the database
- db path: `./ssDNA_tool/DB` (rep and caps)

---

## 1: Dereplication (Clustering)

If you want to cluster (dereplicate) the metagenomic sequences, run:

```bash
python ./ssDNA_tool/ssDNA_annotator/modules/cluster.py \
     -i /path/to/my_sequence.fa \
     -o ./output_clusters
```

### Expected Outputs
The output directory (`output_clusters/`) will contain:
```pgsql
output_clusters/.
├── clusters.tsv           # File with the first column as the cluster representative
├── cluster_sequences.fa     # Clustered FASTA file
├── ani_results.tsv          # Pairwise ANI (average identity) results
├── blast_results.tsv        # BLAST (n/p) results
└── clustering.log           # Module log file
```

## 2: Build a Phylogenetic Tree
### a) Sequence Alignment
You have three options for performing the alignment:

##### Option 1: Use Only Your Sequences
```bash
python ./ssDNA_tool/ssDNA_annotator/modules/align.py  \
                      --threads 24 \
                      --input_fasta /path/to/my_sequence.fa \
                      -d path/to/output/directory
```

##### Option 2: Option 2: Use the Provided Database
Select sequences by family (e.g., `Circoviridae`, `Microviridae`, etc.) or use `all` to include all families.

If using the tool's DB the `path/to/database` should be:
- `./ssDNA_tool/tree/main/DB/caps` for CAPs
- `./ssDNA_tool/tree/main/DB/reps` for REPS

```bash
python ./ssDNA_tool/ssDNA_annotator/modules/align.py  --threads 24 \
                    --input_fasta /path/to/my_sequence.fa \
                    --db_family [Circoviridae Microviridae ...] or [all] \
                    --db_path path/to/database \
                    -d path/to/output/directory
```
##### Option 3: Use Your Own Database
Change the database directory and specify the database file name with `--db_family`.
```bash
python ./ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
       --input_fasta /path/to/my_sequence.fa \
       --db_family name_of_my_db.fa \
       --db_path path/to/my/own/database \
       -d path/to/output/directory
```

##### Alignment Outputs
The alignment output directory (`output/`) will include:
```pgsql
output/.
├── alignment.log                            # Log file
├── metadata.csv                             # Metadata with protein_id, protein_description, family, scientific_name, protein_name, and source.
├── my_sequences_aligned_sequences.fasta     # Sequence alignment FASTA file
├── my_sequences_aligned_trimmed_sequences.fasta  # Trimmed sequence alignment FASTA file
└── my_sequences_merged.fasta                # Merged FASTA file (all sequences used in alignment and trimming)
```

### b) Phylogenetic Tree Construction
Use the (trimmed) alignment file to build the phylogenetic tree:

```bash
python ./ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
       -i ./output/my_sequences_aligned_trimmed_sequences.fasta \
       -d path/to/output/directory/tree

```

#### Tree Building Outputs
The tree output directory (`output/tree/`) will contain:
```pgsql
output/tree/.
├── my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv  
│    # Table with modified sequence names (spaces replaced by "_" to avoid issues in downstream analyses)
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta  
│    # Sequences used to build the phylogenetic tree
├── build_tree.log                           # Module log file
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.bionj  
│    # BIONJ format tree (an improved version of the NJ algorithm)
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.ckp.gz  
│    # IQ-tree checkpoint file
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.contree  
│    # Consensus tree with branch supports (branch lengths optimized on the original alignment)
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.log  
│    # Log file for the entire tree-building run
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.mldist  
│    # Pairwise distance matrix (for downstream analyses)
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.model.gz  
│    # IQ-TREE model checkpoint file
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.splits.nex  
│    # Support values (percentages) for splits, computed from bootstrap trees
├── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.iqtree  
│    # Main IQ-tree report file (self-readable; see computational results)
└── my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.treefile  
     # ML tree in NEWICK format (compatible with viewers like FigTree or iTOL)
```

### c) Plotting the Tree
A plotting module (using [ggtree](https://yulab-smu.top/treedata-book/)) is available to visualize the tree. Note that the script provides automated plotting with some limitations regarding annotation and coloring. For more detailed tree customization, consider using tools like FigTree or iTOL.

 - layout can be one of `'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'`
 - branch_length: if `none` draw cladogram
 - open_angle only for `fan` layout
 - offset: tiplab offset, horizontal adjustment to nudge tip labels, defaults to 0.14
 - tip_label: name of the color group (default is `family` which is based on the [build_tree](https://github.com/ricrocha82/ssDNA_tool/blob/main/ssDNA_annotator/modules/build_tree.py) module metadata)
 - fig_width and fig_height are based on [ggsave](https://ggplot2.tidyverse.org/reference/ggsave.html) function in R

```bash
# Use the teefile from IQ-TREE
python ./ssDNA_tool/ssDNA_annotator/modules/plot_tree.py \
			--tree my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.treefile \
			--outdir ./output/tree \
			--metadata_1 ./output/metadata.csv \
			--metadata_2 ./output/tree/my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv \
			--layout rectangular --branch_length branch.length \
			--open_angle 0 --offset 0.15 \
			--tip_label family \
			--fig_width 20 --fig_height 15 \
			--plot_name my_custom_tree.pdf

# or use the distance table from IQ-TREE
python ./ssDNA_tool/ssDNA_annotator/modules/plot_tree.py \
			--dist_matrix my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.mldist \
			--outdir /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree \
			--metadata_1 ./output/metadata.csv \
			--metadata_2 ./output/tree/my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv \
			--layout rectangular --branch_length branch.length \
			--open_angle 0 --offset 0.15 \
			--tip_label family \
			--fig_width 20 --fig_height 15 \
			--plot_name my_custom_tree_dist.pdf
```
## 3: MOTIF
This module lets you search for a motif in a FASTA file (using the regex pattern provided) and optionally split the sequences at the motif occurrence. It logs every step to a log file named `motif.log` in the specified output directory.

The additional --split flag instructs the module to:
- Read the generated motif positions.
- Split each sequence into two parts:
     - The part before the motif.
     - The part from the motif onward.
- Write the split sequences to two FASTA files (named with a _1.fasta and _2.fasta suffix) in the specified output directory.

```bash
# Running Without Splitting
python ./ssDNA_tool/ssDNA_annotator/modules/motif.py \
    -i my_sequences.fasta \
    -d ./output/motif \
    -p "[GA].{4}GK[TS]"

# Running With Splitting
python ./ssDNA_tool/ssDNA_annotator/modules/motif.py \
    -i my_sequences.fasta \
    -d ./motif_split \
    -p "[GA].{4}GK[TS]" \
    --split

# if using an aligment file, you may want to remove gaps (-) before detecting patterns
python ./ssDNA_tool/ssDNA_annotator/modules/motif.py \
    -i my_sequences_aligned_trimmed_sequences.fasta \
    -d ./motif_split \
    -p "[GA].{4}GK[TS]" \
    --split --remove-gaps
```

## 4: make Sequence log
This module generates sequence logos from FASTA files or motif detection tables. It supports splitting the figure by a metadata column and automatically detects protein vs. nucleotide sequences. A log file (`seq_logo.log`) is automatically created in the output directory.
```bash
# Basic Sequence Logo Generation using the output of the MOTIF module
python seq_logo.py -tb pattern_positions.txt -o output_dir --output_name logo.pdf

# Generating Sequence Logo from a FASTA File
python seq_logo.py -f sequences.fasta -o output_dir --output_name fasta_logo.pdf
```
- -tb pattern_positions.txt: Input motif table (e.g., from seqkit locate).
- -o output_dir: Output directory for the generated logo.
- --output_name logo.pdf: Name of the output sequence logo.

```bash
# Splitting the Figure by Group Labels (Metadata)
python seq_logo.py -tb pattern_positions.txt -o output_dir \
    --output_name logo_split.pdf --split \
    --metadata metadata.csv --ncol 2 --group_label family
```
- --split: Enables splitting based on a metadata column.
- --metadata metadata.csv: Metadata file (must contain seqID and group labels).
- --ncol 2: Number of columns in the split figure.
- --group_label family: Column name in metadata.csv for grouping.

example of the metadata
```pgsql
seqID,family
seq1,Family_A
seq2,Family_B
seq3,Family_A
```