# Motif Discovery

CRESSENT provides comprehensive motif discovery capabilities through two complementary approaches: regex-based pattern matching for known motifs and de novo motif discovery using MEME. Additionally, it integrates with the Prosite database for functional annotation of discovered motifs.

```{image} _static/figures/fig_module_phylogenetic.png
:width: 800
:class: no-scaled-link
:align: center
```

## Overview

The motif discovery module combines:
- **Pattern-based searching** using regex patterns with seqkit
- **De novo motif discovery** using MEME for finding unknown patterns
- **Functional annotation** via ScanProsite for protein sequences
- **Visualization** through sequence logos and motif maps

## Regex-based Pattern Matching

### Basic Usage

Search for specific motif patterns in your sequences:

```bash
cressent motif -i sequences.fasta -p "TAGTATTAC" -o output_dir
```

### Key Features

- **Flexible pattern matching** using regex syntax
- **Gap handling** with optional gap removal
- **Position tracking** with detailed coordinate information
- **Sequence splitting** at motif positions

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-p, --pattern` | Sequence pattern (regex) for motif searching | Required |
| `-n, --table_name` | Name of output table file | `pattern_positions.txt` |
| `--remove-gaps` | Remove gaps before searching | False |
| `--split-sequences` | Split sequences at motif positions | False |

### Example: CRESS Virus Nonanucleotide

```bash
# Search for nonanucleotide motif in CRESS viruses
cressent motif -i cress_genomes.fasta \
               -p "TAGTATTAC" \
               -n nona_positions.txt \
               -o motif_analysis/
```

### Output Files

- **Position table**: Tab-delimited file with motif locations
- **Split sequences**: Optional FASTA files split at motif positions
- **Log file**: Detailed analysis log

## De Novo Motif Discovery

### MEME Integration

Discover unknown motifs using the MEME suite:

```bash
cressent motif_discovery -i sequences.fasta \
                        -o meme_output/ \
                        -nmotifs 3 \
                        -minw 6 \
                        -maxw 12
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-nmotifs` | Number of motifs to find | 1 |
| `-minw` | Minimum motif width | 5 |
| `-maxw` | Maximum motif width | 10 |
| `--meme_extra` | Additional MEME arguments | None |
| `--scanprosite` | Run ScanProsite analysis | False |

### MEME Output Processing

The module automatically processes MEME results to generate:

1. **Consensus table** (`consensus_table.csv`)
2. **Detailed motif table** (`motif_table.csv`) with:
   - Sequence IDs and matched regions
   - Motif positions and orientations
   - Regular expressions for each motif
3. **EPS visualization files** (organized in `eps_files/`)

### Example Output Structure

```
motif_discovery_output/
├── meme.html                 # MEME results webpage
├── meme.xml                  # Machine-readable results
├── consensus_table.csv       # Motif consensus sequences
├── motif_table.csv          # Detailed motif matches
└── eps_files/               # MEME visualization files
    ├── logo1.eps
    └── logo2.eps
```

## ScanProsite Integration

### Functional Annotation

For protein sequences, automatically annotate motifs with known functions:

```bash
cressent motif_discovery -i proteins.fasta \
                        -o output/ \
                        --scanprosite
```

### Features

- **Automatic sequence type detection** (DNA vs protein)
- **Prosite database querying** for functional annotations
- **Rate limiting** to respect server resources
- **Comprehensive results** with functional descriptions

### ScanProsite Output

- **Results table** (`scanprosite_results.csv`) with:
  - Signature accessions and descriptions
  - Pattern positions and scores
  - Functional annotations

## Combined Workflow

### Complete Motif Analysis

```bash
# Step 1: Search for known patterns
cressent motif -i sequences.fasta \
               -p "TAGTATTAC" \
               --generate-logo \
               -o analysis/

# Step 2: Discover novel motifs
cressent motif_discovery -i sequences.fasta \
                        -o analysis/meme/ \
                        -nmotifs 5 \
                        -minw 8 \
                        -maxw 15

# Step 3: Functional annotation (for proteins)
cressent motif_discovery -i proteins.fasta \
                        -o analysis/prosite/ \
                        --scanprosite
```

## Visualization Integration

### Sequence Logos

Generate publication-ready sequence logos:

```bash
cressent motif -i sequences.fasta \
               -p "MOTIF_PATTERN" \
               --generate-logo \
               --plot-title "CRESS Nonanucleotide" \
               --width 12 \
               --height 8 \
               -o output/
```

### Motif Mapping

Create genome-wide motif distribution maps:

```bash
cressent motif_map_viz -f motif_table.csv \
                      -o visualization/ \
                      --format auto
```

## Advanced Features

### Custom MEME Parameters

```bash
cressent motif_discovery -i sequences.fasta \
                        -o output/ \
                        --meme_extra -mod "zoops" -revcomp -dna
```

### Grouped Analysis

Analyze motifs by sequence groups:

```bash
cressent motif -i sequences.fasta \
               -p "PATTERN" \
               --generate-logo \
               --split-logo \
               --metadata groups.csv \
               --group-label "virus_family" \
               --ncol 2 \
               -o grouped_analysis/
```

## Best Practices

### Pattern Design

1. **Use IUPAC codes** for ambiguous positions
2. **Test patterns** on known sequences first
3. **Consider reverse complements** for DNA sequences

### De Novo Discovery

1. **Optimize motif number** based on sequence complexity
2. **Adjust width parameters** for expected motif sizes
3. **Use appropriate background models** for your sequence type

### Performance Tips

1. **Remove gaps** if not biologically relevant
2. **Filter sequences** by length or quality
3. **Use smaller datasets** for initial parameter testing

## Troubleshooting

### Common Issues

**No motifs found:**
- Check pattern syntax
- Verify sequence format
- Consider case sensitivity

**MEME fails:**
- Ensure sequences are aligned if needed
- Check minimum sequence requirements
- Verify MEME installation

**ScanProsite timeout:**
- Reduce sequence number
- Check internet connection
- Retry with rate limiting

### Error Messages

The module provides detailed logging to help diagnose issues:

```
INFO - Processing sequence example_seq: nucleotide sequence (length: 2847)
WARNING - Skipping sequence example_seq2: motif not found
ERROR - MEME command failed: insufficient sequences
```