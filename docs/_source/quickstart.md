# Quickstart

If you are running CRESSENT for the first time, follow the steps below to learn how to perform a complete ssDNA virus analysis. This guide uses Naryaviridae family sequences as an example from [Dai et al](https://www.researchsquare.com/article/rs-6208723/v1), following the workflow demonstrated in our test dataset.

## Setting up the environment

First, activate the CRESSENT conda environment:

```bash
conda activate cressent
```

## Preparing your data

For this tutorial, organize your data directory as follows:

```
data/
├── nary_genomes.fasta       # Nucleotide genome sequences
├── nary_caps.faa            # Capsid protein sequences  
├── nary_reps.faa            # Replication protein sequences
└── nary_proteins.faa        # All protein sequences
```

Since this example uses only 3 sequences, we will skip preprocessing steps like clustering and contamination detection that are typically used for larger datasets.

## Step 1: Genome Analysis

### Nucleotide Sequence Alignment

First, align the nucleotide sequences for recombination analysis:

```bash
cressent align \
    --threads 24 \
    --input_fasta data/nary_genomes.fasta \
    -o output/genome_align
```

### Recombination Detection

Detect recombination events using all available methods:

```bash
cressent recombination \
    -i output/genome_align/nary_genomes_aligned_trimmed_sequences.fasta \
    -o output/recombination \
    -f recomb_results.csv \
    --all
```

### Analyzing Recombination Results

```python
import pandas as pd

# Read recombination results
df = pd.read_csv("output/recombination/recomb_results.csv")

# Filter significant events (p-value < 0.05)
significant = df[df['Pvalue'] < 0.05]

# Count methods detecting each recombinant
method_counts = significant.groupby('Recombinant')['Method'].value_counts().unstack(fill_value=0)
method_counts['Total'] = method_counts.sum(axis=1)

print("Recombination events detected by multiple methods:")
print(method_counts[method_counts['Total'] >= 3])
```

Expected output showing recombination events detected across all genomes:

| Recombinant | 3Seq | Bootscan | Chimaera | MaxChi | RDP | Total |
|-------------|------|----------|----------|--------|-----|-------|
| Genome_1    | 0    | 2        | 1        | 1      | 2   | 6     |
| Genome_2    | 0    | 3        | 1        | 4      | 9   | 17    |
| Genome_3    | 1    | 1        | 0        | 0      | 1   | 3     |

## Step 2: Capsid Protein Analysis

### Sequence Clustering and Alignment

For larger datasets, you would typically cluster sequences first:

```bash
# Skip clustering for small datasets, but clean sequence names
cressent cluster \
    -i data/nary_caps.faa \
    -o output/caps/cluster
```

### Database-Integrated Alignment

Align capsid sequences with the Naryaviridae database for phylogenetic context:

```bash
cressent align \
    --threads 24 \
    --input_fasta data/nary_caps.faa \
    --db_family "Naryaviridae" \
    --protein_type caps \
    --db_path /path/to/databases \
    -o output/caps/align_family
```

### Phylogenetic Tree Construction

Build a phylogenetic tree using an appropriate evolutionary model:

```bash
cressent build_tree \
    -i output/caps/align_family/nary_caps_aligned_trimmed_sequences.fasta \
    -o output/caps/tree \
    -m Q.pfam+F+G4
```

### Tree Visualization

Create publication-ready tree visualizations:

```bash
# Basic circular tree
cressent plot_tree \
    --tree output/caps/tree/nary_caps_aligned_trimmed_sequences.treefile \
    -o output/caps/tree \
    --metadata_1 output/caps/align_family/metadata.csv \
    --metadata_2 output/caps/tree/nary_caps_aligned_trimmed_sequences_sanitized_name_table.tsv \
    --layout circular \
    --offset 0.15 \
    --fig_width 20 --fig_height 15 \
    --plot_name nary_caps_tree.pdf

# Tree with alignment visualization
cressent plot_tree \
    --tree output/caps/tree/nary_caps_aligned_trimmed_sequences.treefile \
    -o output/caps/tree \
    --metadata_1 output/caps/align_family/metadata.csv \
    --metadata_2 output/caps/tree/nary_caps_aligned_trimmed_sequences_sanitized_name_table.tsv \
    --alignment output/caps/align_family/nary_caps_aligned_trimmed_sequences.fasta \
    --layout rectangular \
    --plot_tips False \
    --plot_name nary_caps_tree_with_alignment.pdf
```

### De Novo Motif Discovery

Discover conserved motifs in capsid proteins:

```bash
cressent motif_discovery \
    -i data/nary_caps.faa \
    -o output/caps/motif_discovery \
    -nmotifs 5 -minw 6 -maxw 10 \
    --meme_extra "-mod zoops -evt 0.05" \
    --scanprosite
```

Visualize the discovered motifs:

```bash
cressent motif_map_viz \
    -f output/caps/motif_discovery/scanprosite_results.csv \
    -o output/caps/motif_discovery
```

## Step 3: Replication Protein Analysis

### Alignment and Tree Construction

```bash
# Align replication proteins with database
cressent align \
    --threads 24 \
    --input_fasta data/nary_reps.faa \
    --db_family "Naryaviridae" \
    --protein_type reps \
    --db_path /path/to/databases \
    -o output/reps/align_family

# Build phylogenetic tree
cressent build_tree \
    -i output/reps/align_family/nary_reps_aligned_trimmed_sequences.fasta \
    -o output/reps/tree \
    -m Q.yeast+G4
```

### Known Motif Analysis

Search for the Walker A motif and split sequences at this position:

```bash
cressent motif \
    -i output/reps/align_family/nary_reps_aligned_trimmed_sequences.fasta \
    -o output/reps/motif \
    -p ".{5}GK[TS].{4}" \
    --remove-gaps \
    --split-sequences \
    --generate-logo --split-logo --ncol 2 \
    --metadata output/reps/align_family/metadata.csv \
    --group-label family
```

(understanding-quickstart-outputs)=
## Understanding the outputs

Let's examine the key files generated by our analysis:

### Recombination Results

The `recomb_results.csv` file contains detailed information about detected recombination events:

| Method | Recombinant | Major_Parent | Minor_Parent | Breakpoint_Start | Breakpoint_End | Pvalue |
|--------|-------------|--------------|--------------|------------------|----------------|--------|
| RDP    | Genome_A    | Genome_B     | Genome_C     | 245             | 678            | 0.023  |
| Bootscan | Genome_A  | Genome_B     | Genome_C     | 240             | 685            | 0.034  |

Events detected by multiple methods (≥3) are considered highly reliable.

### Phylogenetic Trees

The `.treefile` contains the phylogenetic tree in Newick format:

```
((Nary_001:0.1234,Nary_002:0.0987):0.0543,(Nary_003:0.2134,Database_seq:0.1876):0.0621);
```

### Metadata Files

The `metadata.csv` file provides context for each sequence:

| protein_id | protein_description | family | scientific_name | source |
|------------|-------------------|--------|-----------------|--------|
| Nary_001 | Rep protein | Naryaviridae | Naryavirus sp. | input |
| DB_seq_1 | Rep protein | Naryaviridae | Reference strain | database |

## Step 4: Advanced Domain Analysis

### Domain-Specific Phylogenetics

After motif splitting, analyze individual protein domains:

```bash
# Align and build tree for domain 1 (Helicase)
cressent align \
    --threads 24 \
    --input_fasta output/reps/motif/split_sequences_1.fasta \
    -o output/reps/domain1

cressent build_tree \
    -i output/reps/domain1/split_sequences_1_aligned_trimmed_sequences.fasta \
    -o output/reps/domain1

# Align and build tree for domain 2 (Endonuclease)
cressent align \
    --threads 24 \
    --input_fasta output/reps/motif/split_sequences_2.fasta \
    -o output/reps/domain2

cressent build_tree \
    -i output/reps/domain2/split_sequences_2_aligned_trimmed_sequences.fasta \
    -o output/reps/domain2
```

### Tanglegram Comparison

Compare phylogenies between different protein domains:

```bash
cressent tanglegram \
    --tree1 output/reps/domain1/split_sequences_1_aligned_trimmed_sequences.treefile \
    --tree2 output/reps/domain2/split_sequences_2_aligned_trimmed_sequences.treefile \
    --label1 "Helicase Domain" \
    --label2 "Endonuclease Domain" \
    --output output/reps/comparison \
    --name_tanglegram "domain_comparison.pdf" \
    --height 11 --width 30
```

### Sequence Logo Generation

Create sequence logos for conserved functional domains:

```bash
# Helicase domain logo
cressent seq_logo \
    -i output/reps/domain1/split_sequences_1_aligned_trimmed_sequences.fasta \
    -o output/reps/domain1 \
    --output_name helicase_logo.pdf \
    --method bits --width 15

# Endonuclease domain logo
cressent seq_logo \
    -i output/reps/domain2/split_sequences_2_aligned_trimmed_sequences.fasta \
    -o output/reps/domain2 \
    --output_name endonuclease_logo.pdf \
    --method prob --width 15
```

## Quality Assessment

### Evaluating Results

**Recombination Analysis:**
- Focus on events detected by ≥3 methods
- Consider p-values < 0.05 as significant
- Manually inspect alignment around breakpoints

**Phylogenetic Analysis:**
- Bootstrap values ≥70% indicate strong support
- Check for reasonable branch lengths
- Ensure biological relevance of groupings

**Motif Analysis:**
- Verify motifs match known functional domains
- Check conservation across sequences
- Validate with literature when possible

## Complete Analysis Script

Here's a complete script that runs the entire analysis:

```bash
#!/bin/bash

# CRESSENT Naryaviridae Analysis Pipeline
echo "Starting CRESSENT analysis of Naryaviridae sequences..."

# Create output directories
mkdir -p output/{genome_align,recombination,caps,reps}

# 1. Genome-level analysis
echo "Step 1: Analyzing genome sequences..."
cressent align --threads 24 --input_fasta data/nary_genomes.fasta -o output/genome_align
cressent recombination -i output/genome_align/nary_genomes_aligned_trimmed_sequences.fasta \
    -o output/recombination -f recomb_results.csv --all

# 2. Capsid protein analysis
echo "Step 2: Analyzing capsid proteins..."
cressent align --threads 24 --input_fasta data/nary_caps.faa \
    --db_family "Naryaviridae" --protein_type caps --db_path databases/ \
    -o output/caps/align_family

cressent build_tree -i output/caps/align_family/nary_caps_aligned_trimmed_sequences.fasta \
    -o output/caps/tree -m Q.pfam+F+G4

cressent plot_tree --tree output/caps/tree/nary_caps_aligned_trimmed_sequences.treefile \
    -o output/caps/tree --layout circular --plot_name caps_tree.pdf

# 3. Replication protein analysis
echo "Step 3: Analyzing replication proteins..."
cressent align --threads 24 --input_fasta data/nary_reps.faa \
    --db_family "Naryaviridae" --protein_type reps --db_path databases/ \
    -o output/reps/align_family

cressent motif -i output/reps/align_family/nary_reps_aligned_trimmed_sequences.fasta \
    -o output/reps/motif -p ".{5}GK[TS].{4}" --split-sequences --generate-logo

# 4. Domain comparison
echo "Step 4: Comparing protein domains..."
cressent align --threads 24 --input_fasta output/reps/motif/split_sequences_1.fasta \
    -o output/reps/domain1
cressent build_tree -i output/reps/domain1/split_sequences_1_aligned_trimmed_sequences.fasta \
    -o output/reps/domain1

cressent align --threads 24 --input_fasta output/reps/motif/split_sequences_2.fasta \
    -o output/reps/domain2
cressent build_tree -i output/reps/domain2/split_sequences_2_aligned_trimmed_sequences.fasta \
    -o output/reps/domain2

cressent tanglegram --tree1 output/reps/domain1/*.treefile \
    --tree2 output/reps/domain2/*.treefile --label1 "Helicase" --label2 "Endonuclease" \
    --output output/reps/comparison --name_tanglegram "domain_comparison.pdf"

echo "Analysis complete! Results are in the 'output' directory."
echo "Key files:"
echo "- Recombination: output/recombination/recomb_results.csv"
echo "- Trees: output/*/tree/*.treefile"
echo "- Visualizations: output/*/tree/*.pdf"
```

## Next Steps

After completing this quickstart tutorial:

1. Explore [individual modules](modules/align.md) for advanced parameter customization
2. Learn about the complete [analysis pipeline](pipeline.md) for complex workflows  
3. Check the [FAQ](faq.md) for common questions and troubleshooting
4. Try analysis with your own ssDNA virus sequences

This quickstart provides a solid foundation for using CRESSENT effectively. The modular design allows you to adapt this workflow for different viral families and research questions.