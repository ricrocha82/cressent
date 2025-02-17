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

---

## Step 1: Dereplication (Clustering)

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

## Step 2: Build a Phylogenetic Tree
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
 - tip_label: name of the color group (default is `family` which is based on the build_tree module metadata)
 - fig_width and fig_height are based on [ggsave](https://ggplot2.tidyverse.org/reference/ggsave.html) function in R

```
# Use the teefile from IQ-TREE
python ./ssDNA_tool/ssDNA_annotator/modules/plot_tree.py \
			--tree=my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.treefile \
			--outdir=./output/tree \
			--metadata_1=./output/metadata.csv \
			--metadata_2=./output/tree/my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv \
			--layout=rectangular --branch_length=branch.length \
			--open_angle=0 --offset=0.15 \
			--tip_label=family \
			--fig_width=20 --fig_height=15 \
			--plot_name=my_custom_tree.pdf

# or use the distance table from IQ-TREE
python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/plot_tree.py \
			--dist_matrix=my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.mldist \
			--outdir=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree \
			--metadata_1=./output/metadata.csv \
			--metadata_2=./output/tree/my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv \
			--layout=rectangular --branch_length=branch.length \
			--open_angle=0 --offset=0.15 \
			--tip_label=family \
			--fig_width=20 --fig_height=15 \
			--plot_name=my_custom_tree_dist.pdf
```
