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
  ()...)

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

If you want to adjust the sequence to begin with determined conserved nonanucleotide sequence run

```bash
python ./ssDNA_tool/ssDNA_annotator/modules/adjust_seq.py  \
                      -i /path/to/my_sequence.fa \
                      -o /path/to/output_directoryu # default is the current directory
                      -m "ATCG" # default: TAGTATTAC
```

If the pattern is found, the output is a fasta file with each sequence beginning with the sequence `-m`

##### Option 2: Use the Database that comes with the tool
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

metadata.csv
```pgsql
protein_id,protein_description,family,scientific_name,protein_name,source
seq1,description_seq1,Family_A,scientific_name1,protein_name1,my_seqs
seq2,description_seq2,Family_B,scientific_name2,protein_name2,db
seq3,description_seq3,Family_A,scientific_name3,protein_name3,my_seqs
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

my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv
```pgsql
Original_Name	Sanitized_Name
Zebra finch circovirus|Circoviridae	Zebra_finch_circovirus_Circoviridae
Bat associated cyclovirus 11|Circoviridae	Bat_associated_cyclovirus_11_Circoviridae
Bat associated cyclovirus 11|Circoviridae	Bat_associated_cyclovirus_11_Circoviridae_1
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

# you can also use the align file to plot the tree with the alginmet
python ./ssDNA_tool/ssDNA_annotator/modules/plot_tree.py \
			--tree my_sequences_aligned_trimmed_sequences_sanitized_sequences.fasta.treefile \
			--outdir ./output/tree \
            --alignment ./output/my_sequences_aligned_trimmed_sequences.fasta \
			--metadata_1 ./output/metadata.csv \
			--metadata_2 ./output/tree/my_sequences_aligned_trimmed_sequences_sanitized_name_table.tsv \
			--layout rectangular --branch_length branch.length \
			--open_angle 0 --offset 0.15 \
			--tip_label family \
			--fig_width 20 --fig_height 15 \
			--plot_name my_custom_tree.pdf 
            
```

### d) Plotting a Tanglegram (or “cophylo plot”) 

```bash
/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/tanglegram.py \
    --tree1 tree1.treefile \
    --tree2 tree2.treefile \
    --label1 my_tree1 \
    --label2 my_tree2 \
    --output .path/to/output/ \
    --name_tanglegram "my_tanglegram.pdf" \
    --height 11 \
    --width 20 \
    --lab_cex 2 # size of lables
```


## 3: MOTIF

### 3.1: Using Regex
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

The tree output directory (`output/motif/`) will contain:
```pgsql
output/motif/.
├── logo.pdf 
├── logo_split.pdf
│  # pdf figures
├── motif.log
├── pattern_positions.txt 
│  # filte containg the pattern positions (seqID patternName pattern strand start end matched)
├── seq_logo.log
│  # log file
│  # sequences if option (--split or/and --remove-gaps) are selected
├── ungapped_sequences.fasta
├── split_sequences_1.fasta
└── split_sequences_2.fasta

```

pattern_positions.txt 
```pgsql
seqID	patternName	pattern	strand	start	end	matched
Rhynchosia golden mosaic Havana virus-[Cuba:Havana:28:2007]|Geminiviridae	[GA].{4}GK[TS]	[GA].{4}GK[TS]	+	235	242	GPSRTGKT
Macroptilium mosaic Puerto Rico virus-[Bean]|Geminiviridae	[GA].{4}GK[TS]	[GA].{4}GK[TS]	+	235	242	GPSRTGKT
Cowpea golden mosaic virus-[Nigeria]|Geminiviridae	[GA].{4}GK[TS]	[GA].{4}GK[TS]	+	235	242	GESRTGKT
(...)
```

### 3.1: _De novo_ Motif Discovery
Pipeline designed to discover __de novo__ motifs from an input FASTA file using [MEME](https://pubmed.ncbi.nlm.nih.gov/7584402/).

It is optional to use [ScanProsite](https://pubmed.ncbi.nlm.nih.gov/16845026/) to scan the protein sequences against the PROSITE collection of motifs

### Basic Options

- `-i`, `--fasta`:The input FASTA file containing the sequences to be analyzed (required). 
- `-o`, `--output`: output directory (optional; default: current directory).
- `-nmotifs`: The number of motifs to discover using MEME. (optional; default: 5).
- `-minw`: The minimum motif width that MEME should conside (optional; default: 6).
- `-maxw`: The maximum motif width that MEME should conside (optional; default: 50).
- `--meme_extra`: Additional MEME arguments provided as key-value pairs. This allows users to pass extra options to MEME without modifying the code (optional).
- `--scanprosite`: Run both MEME and ScanProsite

```bash
# basic
python ./ssDNA_annotator/modules/motif_disc.py <input_fasta_file> -o <output_directory> [-nmotifs N] [-minw MINW] [-maxw MAXW] [--meme_extra KEY VALUE ...] --scanprofite

# only MEME
python ./ssDNA_annotator/modules/motif_disc.py \
                    -i ./seq_meme.fa \
                    -o ./motif_disc \
                    -nmotifs 5 -minw 6 -maxw 50 \
                    --meme_extra -mod zoops -evt 0.05

# MEME + ScanProsite
python ./ssDNA_annotator/modules/motif_disc.py \
                    -i ./seq_meme.fa \
                    -o ./motif_disc \
                    -nmotifs 5 -minw 6 -maxw 50 \
                    --meme_extra -mod zoops -evt 0.05 \
                    --scanprosite

```

#### Outputs
```pgsql
.
├── consensus_table.csv
│  # table with the consensus motif sequences (id,consensus,length,occurrences)
├── logo1.eps
│  # sequence logo in eps format (meme output) for each discovered motif.
├── MEME-1_pwm_matrix.csv
│  # Position Weight Matrix (PWM) matrix for each discovered motif.
├── meme.html # meme output
├── meme.txt # meme output
├── meme.xml # meme output
├── motif_discovery.log
│  # log file
├── motif_table.csv
│  # table with all the motif sequences and attributes (length,motif_name,pvalue,sequence_id,sequence_name,start,end,strand)
└── scanprosite_results.csv (if --scanprofile selected)
    # ScanProsite table (sequence_ac,start,stop,signature_ac,score,level,sequence_id,sequence_name,prosite_ann)
```

You can plot a gene map using the motif_table.csv as input
```bash
python ./ssDNA_annotator/modules/gene_map.py \
                    -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc/motif_table.csv \
                    -o ./motif_disc/gene_map_motif.pdf \
                    --height 10 --width 10 \
                    -t "Motif Distribution in Genes"
```

## 4: make Sequence logo and more
This module generates sequence logos from FASTA files or motif detection tables. It supports splitting the figure by a metadata column and automatically detects protein vs. nucleotide sequences. A log file (`seq_logo.log`) is automatically created in the output directory.

- -tb pattern_positions.txt: Input motif table (e.g., from seqkit locate).
- -o output_dir: Output directory for the generated logo.
- --output_name logo.pdf: Name of the output sequence logo.
- --split: Enables splitting based on a metadata column.
- --metadata: Metadata file (must contain protein_id and group labels). From the aligment module.
- --ncol 2: Number of columns in the split figure.
- --group_label: Column name in metadata.csv for grouping. for example `family`

```bash
# Basic Sequence Logo Generation using the output of the MOTIF module
python seq_logo.py -tb output/motif/pattern_positions.txt -o output/motif/ --output_name logo.pdf

# Generating Sequence Logo from a FASTA File
python seq_logo.py -f my_sequences.fasta -o output/motif/ --output_name fasta_logo.pdf

# Splitting the Figure by Group Labels (Metadata)
python seq_logo.py \
    -tb output/motif/pattern_positions.txt 
    -o output/motif \
    --output_name logo_split.pdf \
    --split \
    --metadata /output/metadata.csv \
    --ncol 2 --group_label family
```

(...) plot the motifs in pile style ([Motifstack](https://bioconductor.org/packages/devel/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html#motifPiles) package)

## 5: Putative Stem Loop and Iterons annotation
### 5.1 Stem Loop Finder

# StemLoop-Finder Module

Customized version of [StemLoop-Finder](https://journals.asm.org/doi/10.1128/mra.00424-21) designed to detect and annotate stem-loop (DNA hairpin) structures with conserved motifs. This module automates the process of finding candidate stem-loop regions using secondary structure prediction. It leverages the `ViennaRNA` library for folding predictions and scores candidate structures based on how closely they match user-defined ideal stem and loop lengths.

## Features

- **Automated Detection:** Searches DNA sequences for potential stem-loop structures near conserved nonanucleotide motifs.
- **Secondary Structure Prediction:** Uses ViennaRNA’s minimum free energy algorithms to predict the dot-bracket structure of candidate regions.
- **Scoring System:** Calculates a score based on deviations from ideal stem and loop lengths. Lower scores indicate candidates closer to the ideal.
- **Flexible Input:** Can use a specific motif (or select a viral family to use a predefined motif) and adjust the number of flanking nucleotides used for folding.
- **Multiple Output Formats:** Produces a GFF file with genome annotations for each predicted stem-loop, and optionally a CSV file with detailed candidate information.
- **(Optional) Logging:** If integrated with logging (similar to the CRUISE module), a log file can be generated to track the analysis process.

## Input Files

- **FASTA File:** Contains the DNA sequences to be analyzed.  
  *(Example: `input.fasta`)*
- **GFF File:** Contains genome annotations, including previously annotated features. This file is used to integrate or overlay new stem-loop annotations.  
  *(Example: `input.gff`)*

## Command-Line Arguments

The StemLoop-Finder module accepts several arguments to customize its behavior:

- `-i`: Path to the input FASTA file.  
- `--gff_in`: Path to the input GFF file.  
- `--out_gff`: Path for the output GFF file containing annotated stem-loops.  
- `--motif`: Conserved nonanucleotide motif to search within candidate regions (allows ambiguous bases). Default:`nantantan`
- `--output_dir`: Path to the output directory
- `--family`: Specify a CRESS DNA virus family (e.g., geminiviridae, circoviridae, etc.) to use a family-specific motif.
- `--idealstemlen` or `-s`: Ideal stem length to score candidate stem-loops. Default:`11`
- `--ideallooplen` or `-l`: Ideal loop length to score candidate stem-loops. Default:`11`
- `--frame` or `-f`: Number of nucleotides flanking the motif (used as the window for folding prediction). Default:`15`
- `--csv_out`: (Optional) Path for an output CSV file with detailed candidate information (e.g., scores, positions, and predicted structures).

## Running the Module

From the command line, run the StemLoop-Finder module as follows:

```bash
python ./ssDNA_annotator/modules/sl_finder.py \
     -i ./test.fasta \
     --gff_in ./test.gff \
     --output_dir ./sl_finder_output \
     --out_gff test_out.gff \
     --csv_out test.csv 
```

Or, if you prefer to specify a motif and/or a family directly:
```bash
python ./ssDNA_annotator/modules/sl_finder.py \
     -i ./test.fasta \
     --gff_in ./test.gff \
     --output_dir ./sl_finder_output \
     --out_gff test_out.gff \ 
     --motif nantantan \
     --family geminiviridae \
     --idealstemlen 11 --ideallooplen 11 --frame 15
```
#### Output files
```pgsql
output/
├── stemloop_finder.log   # Log file capturing the analysis process
├── output.gff            # GFF file with stem-loop annotations
└── output.csv            # (Optional) CSV file with detailed candidate data
```

### 5.1 Cruise
Customized version of [CRUISE](https://journals.asm.org/doi/10.1128/mra.01123-22) (CRiteria-based Uncovering of Iteron SEquences) to search virus genomes for iterons. The module scans around nonanucleotide features and stem-loop structures from an input GFF file (with associated FASTA sequences) to identify and score candidate iteron repeats. Iteron annotations are then output to a new GFF file, and a detailed log is saved for tracking the analysis process.

#### Features
- **Iteron Detection**: Searches for candidate iteron substrings in the genomic sequence.
- **Scoring and Ranking**: Scores iteron candidates based on distance criteria, length, and proximity to stem-loop features.
- **Annotation Output**: Annotates iteron and stem-loop repeat features into a new GFF file.
- **Logging**: Generates a log file (cruise.log) that records key steps and outcomes.
- **Customizable Parameters**: Accepts several command-line arguments to fine-tune the iteron search.

```bash 
python cruise.py --inputFasta examples/test.fasta \
                 --inputGFF examples/out.gff \
                 --outputGFF examples/finaloutput.gff \
                 --outputDir output/ \
                 --verbose
```

#### Input Files
**Input FASTA File**: Contains the complete genome sequences. (Example: `examples/test.fasta`)
**Input GFF File**: Contains genome annotations including nonanucleotide features and stem-loop features. (Example: `examples/out.gff`)

#### Arguments
The module accepts the following arguments:

- `--inputFasta`: Path for input FASTA file with all sequences.
- `--inputGFF`: Path for associated input GFF file.
- `--outputGFF`: Path for the output GFF file (with iteron annotations).
- `--minLength`: Minimum iteron length (in nucleotides). Default: 5
- `--maxLength`: Maximum iteron length (in nucleotides). Default: 12
- `--range`: Number of base pairs around the nonanucleotide to search. Default: 65
- `--rank`: Use ranking system for iteron candidates (True/False). Default: True
- `--numberTopIterons`: Number of iterons returned in rank order. Default: 5
- `--maxScore`: Maximum score allowed for iterons if ranking is not used. Default: 40
- `--wiggle`: Max difference between iteron length and distance allowed before score detriment. Default: 5
- `--goodLength`: Highest favorable iteron length. Default: 11
- `--doStemLoop`: Whether to annotate stem-loop repeats (True/False). Default: True
- `--doKnownIterons`: Whether to annotate known iterons (True/False). Default: True
- `--maxDist`: Maximum allowed distance between iteron occurrences. Default: 20
- `--bestDist`: Optimal maximum distance between iterons. Default: 10
- `--scoreRange`: Score range between output candidates (used for filtering ranked iterons). Default: 50
- `--outputDir`: Directory to save all output files (GFF, log, etc.). Default: "." (current directory)
- `--verbose`: Enable verbose output (prints progress to the console).

#### Output files
```pqsql
output/
├── cruise.log
└── finaloutput.gff
```

## 6: Recombination detection
The recombination module is designed to detect recombination events in nucleotide sequences using multiple detection methods. It integrates a customized version of [OpenRDP](https://github.com/aglucaci/OpenRDP/tree/master) ([Recombination Detection Program](https://academic.oup.com/ve/article/7/1/veaa087/6020281)) to provide a comprehensive suite of recombination detection algorithms.

When run for the first time, recombination.py will if binary executables exist (3Seq and GENECONV), compile from source if not. Also, it will generate a new 500 × 500 × 500 P-value table. If you don't want to genearate a p-value table you can extract from [here](https://github.com/ricrocha82/ssDNA_tool/blob/main/ssDNA_annotator/modules/openrdp/bin/source_code/myPvalueTable.tar.gz)

```bash
python recombination.py -i <input_alignment> -o <output_file> [options]

# Run all methods
python recombination.py -i aligned_sequences.fasta -o results.csv

# Specify output directory
python recombination.py -i aligned_sequences.fasta -o results.csv -d output_dir

# Run specific methods
python recombination.py -i aligned_sequences.fasta -o results.csv -rdp -maxchi -bootscan

# Use custom configuration
python recombination.py -i aligned_sequences.fasta -o results.csv -c my_config.ini -all
```
### Basic Options

- `-i`, `--input`: Input alignment file in FASTA format (required). 
    - **Input Alignment File**: A multiple sequence alignment in FASTA format
    - Must contain at least 3 sequences
    - All sequences must be of the same length (aligned)
    - Should be in standard nucleotide format (A, T, G, C)
- `-o`, `--output`: Output file for results in CSV format (required)
- `-d`, `--outdir`: Output directory for all files (default: current directory)
- `-c`, `--config`: Configuration file in INI format for OpenRDP parameters
    - **Configuration File**: An INI format file with parameters for each method
    - If not provided, default parameters will be used
    - Example configuration files can be found in the `scripts/default_config.ini`


### Method Selection

- `-rdp`: Run RDP method
- `-threeseq`: Run 3Seq method
- `-geneconv`: Run GENECONV method
- `-maxchi`: Run MaxChi method
- `-chimaera`: Run Chimaera method
- `-bootscan`: Run Bootscan method
- `-siscan`: Run Siscan method
- `-all`: Run all methods (default if no method is specified)

#### Methods

The module implements seven recombination detection methods:

1. **RDP**: Recombination Detection Program
   - Uses a sliding window to identify changes in sequence similarity patterns
   - Good at detecting recent recombination events

2. **3Seq**: [3-Sequence Method](https://mol.ax/software/3seq/)
   - Detects recombination by examining triplets of sequences
   - Tests whether a sequence is a recombinant of two parental sequences
   - Robust to rate heterogeneity

3. **GENECONV**: [Gene Conversion Detection](https://www.math.wustl.edu/~sawyer/geneconv/)
   - Detects gene conversion events by identifying unusually similar sequence fragments
   - Particularly useful for detecting older recombination events

4. **MaxChi**: Maximum Chi-Square
   - Uses chi-square statistic to detect breakpoints
   - Looks for significant differences in proportions of variable sites

5. **Chimaera**: Similar to MaxChi
   - Modified version of MaxChi with different statistical approach
   - Often detects breakpoints missed by other methods

6. **Bootscan**: Bootscanning Method
   - Analyzes phylogenetic signal changes along the sequence
   - Identifies regions where evolutionary relationships change

7. **Siscan**: Sister-Scanning Method
   - Uses statistical tests to identify changes in phylogenetic relationships
   - Particularly good at detecting distant recombination events

### Output Control

- `-quiet`: Suppress console output
- `-verbose`: Enable verbose logging

## 7: GC Content Heatmap
Calculates %G+C content using a sliding window approach.
Output is a heatmap with GC content. Rows as sequence names (left) and % of GC (right).

```bash
python ./ssDNA_annotator/modules/gc_patchiness.py \
                                          -i ./test_fast.fa \
                                        --output_dir ./output \
                                          --output_name gc_heatmap.pdf \
                                        --fig_width 10 --fig_height 20 
```



