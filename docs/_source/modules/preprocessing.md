# Preprocessing

CRESSENT provides comprehensive preprocessing capabilities to prepare sequence data for downstream analysis. The preprocessing module includes dereplication, decontamination, and sequence adjustment tools to ensure high-quality input data.

```{image} _static/figures/fig_module_preprocessing.png
:width: 800
:class: no-scaled-link
:align: center
```

## Overview

The preprocessing module offers three key functionalities:
- **Dereplication** using CD-HIT and clustering algorithms
- **Decontamination** via BLAST-based screening against contaminant databases
- **Sequence adjustment** for motif-based sequence standardization

## Dereplication and Clustering

### Purpose

Remove redundant sequences and group similar sequences to reduce computational burden and improve analysis quality.

### Basic Usage

```bash
cressent cluster -i sequences.fasta \
                 -o clustering_output/ \
                 -t 8 \
                 --min_ani 95.0 \
                 --min_tcov 85.0
```

### Algorithm Workflow

1. **Sequence preprocessing** with name sanitization
2. **BLAST database creation** (nucleotide or protein auto-detected)
3. **All-vs-all BLAST** search
4. **ANI calculation** using anicalc
5. **Clustering** with aniclust
6. **Representative selection** and output generation

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --threads` | Number of CPU threads | 1 |
| `--min_ani` | Minimum average nucleotide identity | 95.0 |
| `--min_tcov` | Minimum target coverage | 85.0 |
| `--min_qcov` | Minimum query coverage | 0.0 |
| `--keep_names` | Keep only first word of sequence IDs | False |

### Output Files

- **Clusters table** (`clusters.tsv`): Representative sequences and members
- **Representative sequences** (`cluster_sequences.fa`): FASTA with cluster representatives
- **Name mapping** (`*_name_table.tsv`): Original to sanitized name mapping
- **BLAST results** (`blast_results.tsv`): All-vs-all comparison results
- **ANI results** (`ani_results.tsv`): Average nucleotide identity calculations

### Example Output

```
clustering_output/
├── clusters.tsv
├── cluster_sequences.fa
├── blast_results.tsv
├── ani_results.tsv
├── renamed_sequences.fasta
├── sequences_name_table.tsv
└── clustering.log
```

### Cluster Table Format

```
Representative_Sequence    Sequences
seq_001                   seq_001,seq_045,seq_123
seq_002                   seq_002,seq_067
seq_003                   seq_003
```

## Decontamination

### Purpose

Screen and remove potential contaminant sequences using BLAST against known contaminant databases.

### Basic Usage

```bash
cressent detect_contamination -i input_sequences.fasta \
                             --db contaminant_database.fasta \
                             -o decontamination_output/ \
                             --output-name clean_sequences
```

### Features

- **Automatic sequence type detection** (nucleotide vs protein)
- **Flexible BLAST parameters** for different stringency levels
- **Comprehensive statistics** and contamination reports
- **Cross-platform compatibility** with robust error handling

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--db` | Contaminant database FASTA file | Required |
| `--seq-type` | Sequence type (nucl/prot) | Auto-detect |
| `--evalue` | BLAST E-value threshold | 1e-10 |
| `--identity` | Minimum percent identity | 90.0 |
| `--coverage` | Minimum query coverage | 50.0 |
| `-t, --threads` | Number of CPU threads | 1 |
| `--keep-temp` | Keep temporary BLAST files | False |

### Algorithm Steps

1. **Input validation** and sequence type detection
2. **BLAST database creation** from contaminant sequences
3. **BLAST search** (BLASTN or BLASTP)
4. **Results filtering** by identity and coverage thresholds
5. **Sequence filtering** and clean sequence output
6. **Statistics generation** and reporting

### Output Files

- **Clean sequences** (`clean_sequences.fasta`): Filtered sequences
- **Statistics** (`clean_sequences_stats.txt`): Contamination summary
- **BLAST results** (`clean_sequences_blast.tsv`): Optional detailed results
- **Log file** (`clean_sequences_decontamination.log`): Process log

### Statistics Report

```
Total sequences: 1000
Identified contaminants: 45
Clean sequences: 955
Contamination rate: 4.50%

Contaminant sequence IDs:
seq_001
seq_042
...
```

### Building Contaminant Databases

Create custom contamination databases from accession lists:

```bash
cressent build_contaminant_db --accession-csv accessions.csv \
                             -o database_output/ \
                             --output-name viral_contaminants \
                             --email your.email@domain.com
```

#### Features

- **NCBI integration** for sequence download
- **Batch processing** with rate limiting
- **Protein extraction** from nucleotide records
- **Metadata generation** for traceability

## Sequence Adjustment

### Purpose

Standardize sequence start positions based on conserved motifs, particularly useful for circular genomes.

### Basic Usage

```bash
cressent adjust_seq -i sequences.fasta \
                   -m "TAGTATTAC" \
                   -o adjusted_output/
```

### Features

- **Automatic sequence type detection** (DNA/RNA/protein)
- **Flexible motif patterns** using regex syntax
- **Circular genome support** with sequence rotation
- **Comprehensive logging** and statistics

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-m, --motif` | Motif pattern for adjustment | TAGTATTAC |
| `-o, --output` | Output directory | Current directory |

### Algorithm

1. **Sequence type detection** based on character composition
2. **Motif searching** using regex patterns
3. **Sequence rotation** to start at motif position
4. **Quality control** and validation
5. **Output generation** with adjusted sequences

### Output Files

- **Adjusted sequences** (`*_motif_adj.fa`): Sequences starting at motif
- **Log file** (`adjust_seq.log`): Processing details and statistics

### Processing Summary

```
Processing Summary:
Total sequences processed: 100
Sequences adjusted: 87
Sequences skipped: 13
Nucleotide sequences: 95
Protein sequences: 5
Unknown sequences: 0
```

## Combined Preprocessing Workflow

### Complete Pipeline

```bash
# Step 1: Sequence adjustment
cressent adjust_seq -i raw_sequences.fasta \
                   -m "TAGTATTAC" \
                   -o preprocessing/

# Step 2: Decontamination
cressent detect_contamination -i preprocessing/*_motif_adj.fa \
                             --db viral_contaminants.fasta \
                             -o preprocessing/ \
                             --output-name decontaminated

# Step 3: Dereplication
cressent cluster -i preprocessing/decontaminated.fasta \
                 -o preprocessing/ \
                 --min_ani 95.0 \
                 --min_tcov 85.0 \
                 -t 8
```

### Quality Control Workflow

```bash
# Conservative decontamination
cressent detect_contamination -i sequences.fasta \
                             --db contaminants.fasta \
                             --identity 85.0 \
                             --coverage 70.0 \
                             -o qc_output/

# Strict clustering
cressent cluster -i qc_output/clean_sequences.fasta \
                 --min_ani 98.0 \
                 --min_tcov 90.0 \
                 --min_qcov 10.0 \
                 -o qc_output/
```

## Best Practices

### Sequence Preparation

1. **Quality filtering** before preprocessing
2. **Appropriate motif selection** for sequence adjustment
3. **Comprehensive contaminant databases** for screening

### Parameter Optimization

1. **ANI thresholds** based on taxonomic level of interest
2. **Coverage parameters** depending on sequence quality
3. **Identity cutoffs** appropriate for contamination type

### Computational Considerations

1. **Thread allocation** based on available resources
2. **Memory requirements** for large datasets
3. **Temporary file management** for disk space

## Troubleshooting

### Common Issues

**No sequences adjusted:**
- Verify motif pattern syntax
- Check sequence format and encoding
- Consider motif orientation

**High contamination rates:**
- Review contamination database completeness
- Adjust identity/coverage thresholds
- Validate input sequence quality

**Clustering failures:**
- Check sequence format compatibility
- Verify BLAST installation
- Monitor memory usage

### Performance Optimization

**Large datasets:**
- Increase thread count
- Use sequence splitting for memory management
- Consider preliminary filtering

**Slow BLAST searches:**
- Reduce database size if possible
- Optimize E-value thresholds
- Use faster BLAST variants

### Error Diagnostics

The preprocessing modules provide detailed logging:

```
INFO - Processing FASTA with keep_names=False
INFO - Sanitized 1000 sequences to output.fasta
INFO - Detected nucleotide sequences (blastn will be used)
WARNING - Binary compatibility issue detected: GLIBC version
ERROR - BLAST search failed: insufficient memory
```