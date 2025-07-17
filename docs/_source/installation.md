# Installation

## Installing CRESSENT

You can install CRESSENT using either a general-purpose package manager (`mamba`, `conda`) or by cloning the repository and setting up the environment manually.

::::{tab-set}

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/) is a package manager that handles all your dependencies for you. To install CRESSENT using Mamba, you need to create a new environment and activate it before running the `cressent` command.

```bash
mamba create -n cressent 
mamba install -c bioconda -c conda-forge -c defaults cressent 
conda activate cressent
```
:::

:::{tab-item} Manual Installation
If you prefer to install CRESSENT manually, you can clone the repository and set up the environment:

```bash
git clone https://github.com/ricrocha82/cressent.git
cd cressent
mamba env create -f cressent_env.yaml
conda activate cressent
```
:::

::::

```{admonition} Dependencies
:class: attention
CRESSENT depends on several external tools that must be available in your PATH:
- [IQ-TREE](http://www.iqtree.org/) for phylogenetic analysis
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for multiple sequence alignment
- [TrimAl](http://trimal.cgenomics.org/) for alignment trimming
- [MEME Suite](https://meme-suite.org/) for motif discovery
- [SeqKit](https://bioinf.shenwei.me/seqkit/) for sequence manipulation
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for sequence similarity searches
- [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) for sequence clustering
- [MMseqs2](https://github.com/soedinglab/MMseqs2) for fast sequence searches
- [R](https://www.r-project.org/) with required packages for visualization
```

## Verification

After installation, verify that CRESSENT is working correctly:

```bash
# Test basic functionality
cressent --help

# Test individual modules
cressent align --help
cressent build_tree --help
cressent plot_tree --help
```

You should see help messages for each command without any errors.

## Running CRESSENT using containers

You can also execute CRESSENT using containerization tools, such as Docker and Podman.

::::{tab-set}

:::{tab-item} Docker
```bash
docker pull your-registry/cressent:latest
```
:::

:::{tab-item} Podman
```bash
podman pull your-registry/cressent:latest
```
:::

::::

To start a CRESSENT container you have to mount a folder from the host system into the container with the `-v` argument. The following command mounts the current working directory (`$(pwd)`) under `/app` inside the container:

::::{tab-set}

:::{tab-item} Docker
```bash
docker run -ti --rm -v "$(pwd):/app" your-registry/cressent:latest align --help
docker run -ti --rm -v "$(pwd):/app" your-registry/cressent:latest build_tree input.fasta output/
```
:::

:::{tab-item} Podman
```bash
podman run -u 0 -ti --rm -v "$(pwd):/app" your-registry/cressent:latest align --help
podman run -u 0 -ti --rm -v "$(pwd):/app" your-registry/cressent:latest build_tree input.fasta output/
```
:::

::::

## Database Setup

Some CRESSENT modules require reference databases. If you want use the custom Db from CRESSENT please, set up the databases after installation:

### Download Pre-built Databases

```bash
# Create database directory
mkdir -p DB
cd DB
```
- [CAPS](https://github.com/ricrocha82/cressent/tree/main/DB/caps): Capsid protein sequences grouped by family
- [REPS](https://github.com/ricrocha82/cressent/tree/main/DB/reps): Replication-associated protein sequences grouped by family
- [Contamination](https://github.com/ricrocha82/cressent/blob/main/DB/contaminant/contaminant_db.fasta) DB
- [tree_models.csv](https://github.com/ricrocha82/cressent/blob/main/DB/taxonomy_accession_number.csv): table with the models computed by ModelFinder from IQ-TREE2 for each Rep and Cap DB grouped by family

The Databse is also available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15981951.svg)](https://doi.org/10.5281/zenodo.15981951)

### Build Custom Databases

Alternatively, build databases for specific [viral families](https://github.com/ricrocha82/cressent/blob/main/DB/taxonomy_accession_number.csv):

```bash
# Build database for specific viral families
cressent db_builder \
    -t taxonomy_accession_number.csv \
    -l Family \
    -s "Circoviridae" "Genomoviridae" \
    -o custom_database \
    -e your.email@example.com
```

## Troubleshooting

### Common Installation Issues

**Conda/Mamba Environment Conflicts**
```bash
# Clean conda cache and try again
conda clean --all
mamba create -n cressent_clean -c conda-forge -c bioconda [packages]
```

**Missing Dependencies**
```bash
# Check if all tools are in PATH
which iqtree
which mafft
which trimal
which meme

# Install missing tools individually
conda install -c bioconda iqtree
```

**R Package Issues**
```bash
# Install R packages manually
R
> install.packages(c("ggplot2", "ape", "dendextend"))
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install(c("ggtree", "Biostrings"))
```

**Permission Errors**
```bash
# Ensure proper permissions for installation directory
chmod +x /path/to/cressent/bin/*
```

### Performance Optimization

**Memory Usage**
- For large datasets, ensure at least 16GB RAM
- Use `--threads` parameter appropriately (typically number of CPU cores)
- Monitor memory usage with `htop` or `top`

### Getting Help

If you encounter installation issues:

1. Check the [FAQ](faq.md) for common solutions
2. Search existing [GitHub issues](https://github.com/ricrocha82/cressent/issues)
3. Create a new issue with:
   - Your operating system and version
   - Complete error messages
   - Installation method used
   - Output of `conda list` or `pip list`
