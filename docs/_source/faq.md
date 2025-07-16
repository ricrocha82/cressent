# Frequently Asked Questions

Here you will find answers to common questions about CRESSENT. If you can't find an answer to your problem here, please open an issue in the [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/ricrocha82/cressent/).

## General Usage

### How can I speed up CRESSENT analysis?

If you want to speed up CRESSENT execution, several options are available:

- **Use sequence clustering**: Run the [cluster](preprocessing.md) module before phylogenetic analysis to remove redundant sequences
- **Optimize thread usage**: Use the `--threads` parameter with appropriate values for your system
- **Disable time-intensive modules**: Skip [recombination](recombination.md) analysis if not needed for your research
- **Use simpler models**: Choose faster evolutionary models (e.g., `-m GTR+G4` instead of `-m MFP`) for quick analysis

Please note that some optimizations may reduce analysis depth or accuracy.

### What input formats does CRESSENT accept?

CRESSENT accepts the following input formats:

- **FASTA files**: Nucleotide and protein sequences (required)
- **Compressed files**: `.gz`, `.bz2`, and `.xz` compression supported
- **GFF files**: For structural analysis modules like [stem loops and iterons](secondary_structure.md)
- **CSV files**: For database construction and metadata integration

All sequence files must be in valid FASTA format with unique sequence identifiers.

### How do I interpret phylogenetic tree support values?

Bootstrap support values in CRESSENT trees indicate confidence in tree topology:

- **≥95%**: Very strong support, high confidence in the relationship
- **70-94%**: Strong support, generally reliable
- **50-69%**: Moderate support, interpret with caution  
- **<50%**: Poor support, relationship is uncertain

For critical evolutionary inferences, focus on relationships with ≥70% bootstrap support.

## Module-Specific Questions

### Why is my alignment empty or of poor quality?

Poor alignment quality can result from several issues:

**Input sequence problems:**
  - Sequences are not homologous
  - Poor sequence quality with many ambiguous nucleotides
  - Sequences of very different lengths

**Solutions:**
  - Verify that input sequences represent the same gene/protein
  - Remove sequences with >10% ambiguous characters
  - Use [cluster](preprocessing.md) to remove highly divergent sequences
  - Adjust `--mafft_ep` parameter for better alignment sensitivity

### Why does tree construction fail with "insufficient data" errors?

Tree construction can fail when:

- **Too few informative sites**: Sequences are too similar or too short
- **Identical sequences**: Remove duplicates using [cluster](preprocessing.md)
- **Poor alignment**: Insufficient overlap between sequences after trimming

**Solutions:**
- Check alignment quality and length
- Use longer sequences or more divergent taxa
- Adjust trimming parameters (`--gap_threshold`)
- Try simpler evolutionary models

### How can I improve recombination detection sensitivity?

To enhance recombination detection:

1. **Use high-quality alignments**: Ensure proper sequence alignment
2. **Include diverse sequences**: Adequate evolutionary distance improves detection
3. **Run all methods**: Use `--all` flag for comprehensive analysis
4. **Verify manually**: Inspect alignment around detected breakpoints
5. **Apply statistical filters**: Focus on events with p-values <0.05

### Why are my sequence logos poorly resolved?

Poor sequence logo quality often indicates:

- **Insufficient sequence conservation**: Not enough conserved positions
- **Poor motif definition**: Pattern may be too broad or specific
- **Alignment gaps**: Excessive gaps in the motif region

**Solutions:**
- Increase sequence number in the analysis
- Refine motif search patterns
- Use `--remove-gaps` option in [motif](motif.md)
- Check alignment quality in the motif region

## Technical Issues

### Why am I getting memory errors?

Memory errors typically occur with:

- Large datasets (>1000 sequences)
- Long alignments (>10,000 positions)
- Database integration with large reference sets

**Solutions:**
- **Reduce dataset size**: Use [cluster](preprocessing.md) for sequence reduction
- **Decrease thread count**: Lower `--threads` parameter to reduce memory usage
- **Increase system RAM**: Consider upgrading hardware for large analyses
- **Process in batches**: Split large datasets into smaller chunks

### Why is the execution stuck during tree building?

Tree construction can be slow due to:

- **Large datasets**: Many sequences or long alignments
- **Complex models**: Model selection (`-m MFP`) is computationally intensive
- **Insufficient CPU**: Single-threaded execution on multi-core systems

**Solutions:**
- Use `--threads` parameter to utilize multiple cores
- Choose specific models (`-m GTR+G4`) instead of model selection
- Monitor system resources to identify bottlenecks
- Consider using a high-performance computing cluster

### Why do binary compilation errors occur?

Some modules compile required binaries automatically:

- **System compatibility**: Binaries may not work on all systems
- **Missing dependencies**: Required compilers or libraries missing
- **Permission issues**: Insufficient write permissions

**Solutions:**
- Ensure `gcc` and `make` are installed
- Check system compatibility (Linux/Unix preferred)
- Run with appropriate permissions
- Install missing system dependencies

## Database and Integration

### Why does database integration fail?

Database integration issues often result from:

- **Missing database files**: Verify database paths and file existence
- **Incompatible family names**: Check spelling and available families
- **Corrupted downloads**: Re-download database files if corrupted

**Solutions:**
- Verify `--db_path` points to correct directory
- Check available families in database documentation
- Re-run database construction if files are corrupted

### How do I handle large-scale comparative analyses?

For large comparative studies:

1. **Use hierarchical clustering**: Reduce dataset complexity first
2. **Employ database integration**: Add reference sequences for context
3. **Optimize computational resources**: Use high-performance computing
4. **Process by viral family**: Analyze related viruses separately
5. **Validate key findings**: Focus detailed analysis on important relationships

## Results Interpretation

### How reliable are the recombination predictions?

Recombination reliability depends on:

- **Multiple method agreement**: Events detected by ≥3 methods are most reliable
- **Statistical significance**: p-values <0.01 indicate strong evidence
- **Biological plausibility**: Results should make biological sense
- **Manual validation**: Visual inspection of alignments confirms events

### What do tanglegram comparisons reveal?

Tanglegrams show congruence between phylogenetic trees:

- **Parallel lines**: Congruent phylogenies, similar evolutionary history
- **Crossed lines**: Incongruent relationships, possible recombination
- **RF scores**: Robinson-Foulds distance quantifies tree differences

High congruence suggests vertical inheritance; incongruence may indicate recombination or different evolutionary pressures.

### How do I validate phylogenetic results?

Validate phylogenetic analyses through:

1. **Bootstrap assessment**: Focus on well-supported relationships (≥70%)
2. **Model adequacy**: Check model fit statistics in log files  
3. **Biological consistency**: Ensure results match known biology
4. **Alternative methods**: Compare with other phylogenetic approaches
5. **Literature comparison**: Verify consistency with published studies

## Best Practices

### What is the recommended workflow for new users?

For new users, we recommend:

1. **Start with tutorial data**: Practice with provided example datasets
2. **Begin with simple analyses**: Single-gene phylogenies before complex workflows
3. **Validate each step**: Check intermediate outputs for quality
4. **Use default parameters**: Optimize parameters only after understanding basics
5. **Consult documentation**: Read module-specific documentation thoroughly

### How should I cite CRESSENT in publications?

Please cite CRESSENT as:

> [Authors]. CRESSENT: A comprehensive toolkit for CRESS DNA virus analysis. *[Journal]* [Year]. DOI: [DOI].

Additionally, cite specific tools used within CRESSENT modules (MAFFT, IQ-TREE, etc.) as appropriate for your analysis.

### What are common analysis pitfalls to avoid?

Avoid these common mistakes:

1. **Skipping quality control**: Always validate input data quality
2. **Ignoring statistical support**: Don't over-interpret poorly supported results
3. **Using inappropriate models**: Choose evolutionary models suitable for your data
4. **Mixing sequence types**: Don't combine nucleotide and protein sequences inappropriately
5. **Inadequate taxon sampling**: Include sufficient diversity for robust analysis

### How do I optimize CRESSENT for different viral families?

Different viral families may require specific approaches:

**Small genomes (e.g., Circoviridae)**:
- Use complete genome sequences when possible
- Consider codon-based alignments for coding sequences
- Pay attention to overlapping genes

**Large genomes (e.g., some ssDNA viruses)**:
- Focus on conserved gene regions
- Use gene-specific databases when available
- Consider computational resource requirements

**Highly divergent families**:
- Use protein sequences instead of nucleotides
- Employ more sensitive alignment parameters
- Consider domain-specific analysis

If you have additional questions not covered here, please check the individual module documentation or contact us through the GitHub repository.