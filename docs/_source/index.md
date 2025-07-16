---
hide-toc: true
---

# CRESSENT

CRESSENT (CRESS DNA Virus Analysis Tool) is a comprehensive bioinformatics pipeline designed for the analysis of ssDNA viruses. It provides state-of-the-art tools for phylogenetic analysis, recombination detection, motif discovery, and functional annotation of CRESS DNA viruses.

```{image} _static/figures/fig_cressent_new.png
:width: 400
:class: no-scaled-link
:align: center
```

::::{grid}
:gutter: 2

:::{grid-item-card}
:columns: 12 12 4 4
**Comprehensive Analysis**
^^^
CRESSENT integrates multiple analysis modules from sequence preprocessing to advanced phylogenetic comparisons.
:::

:::{grid-item-card}
:columns: 12 12 4 4
**Phylogenetic Tools**
^^^
Build publication-ready phylogenetic trees with integrated alignment visualization and domain-specific analysis.
:::

:::{grid-item-card}
:columns: 12 12 4 4
**Motif Discovery**
^^^
Discover and visualize conserved motifs using both known patterns and de novo discovery methods.
:::

::::

## {octicon}`rocket;0.85em` Get started

To start using CRESSENT, read the installation and quickstart guides below. To learn more about specific analysis modules, visit the module pages listed in the sidebar.

:::{card} Installation
:link: installation
:link-type: doc
Instructions on how to install CRESSENT and its dependencies.
:::

:::{card} Quickstart
:link: quickstart
:link-type: doc
Learn how to run CRESSENT and interpret its results with a step-by-step tutorial.
:::

## {octicon}`tools;0.85em` Analysis Modules

CRESSENT provides a modular analysis framework that can be customized for different research needs:

:::{card} Preprocessing
:link: preprocessing
:link-type: doc
Quality control, dereplication, decontamination, and sequence adjustment tools.
:::

:::{card} Alignment
:link: phylogenetic_analysis
:link-type: doc
A module to aligmen your sequences with a custom database.
:::

:::{card} Phylogenetic Analysis
:link: phylogenetic_analysis
:link-type: doc
Sequence alignment, tree building, and comparative phylogenetic analysis.
:::

:::{card} Motif Discovery
:link: motif
:link-type: doc
Pattern-based searching and de novo motif discovery with functional annotation.
:::

:::{card} Secondary Structure Detection
:link: secondary_structure
:link-type: doc
Stem-loop and iteron identification for viral replication elements.
:::

:::{card} Recombination Detection
:link: recombination
:link-type: doc
Comprehensive recombination analysis using multiple detection methods.
:::

:::{card} Visualization Tools
:link: viz
:link-type: doc
Publication-ready figures including sequence logos, trees, and motif maps.
:::

## {octicon}`workflow;0.85em` Pipeline Overview

CRESSENT workflows are designed to be modular and flexible, allowing researchers to combine different analysis modules based on their specific research questions:

### Basic Workflow

1. **Data Preprocessing** - Quality control and sequence preparation
2. **Motif Analysis** - Pattern discovery and functional annotation
3. **Structural Analysis** - Secondary structure detection
4. **Phylogenetic Analysis** - Evolutionary relationships
5. **Visualization** - Publication-ready figures

### Advanced Workflows

- **Comparative Genomics** - Multi-genome analysis with recombination detection
- **Evolutionary Analysis** - Deep phylogenetic analysis with tanglegrams
- **Functional Annotation** - Comprehensive motif and domain analysis

## {octicon}`graph;0.85em` Key Features

### Comprehensive Analysis Pipeline
- **Multi-format support** for FASTA, GFF, and various annotation formats
- **Scalable processing** from single genomes to large datasets
- **Quality control** with contamination detection and sequence validation
- **Reproducible workflows** with detailed logging and parameter tracking

### Advanced Motif Analysis
- **Pattern-based search** using regex with seqkit integration
- **De novo discovery** using MEME for unknown motifs
- **Functional annotation** via ScanProsite database queries
- **Visualization** through information-rich sequence logos

### Structural Biology Tools
- **Stem-loop detection** using ViennaRNA folding algorithms
- **Iteron identification** with CRUISE for replication origins
- **Family-specific analysis** for major CRESS virus groups
- **Validation scoring** based on structural and sequence features

### Phylogenetic Capabilities
- **Multiple alignment** with MAFFT and trimming with TrimAl
- **Tree building** using IQ-TREE with ModelFinder
- **Advanced visualization** with metadata integration
- **Comparative analysis** through tanglegrams and distance matrices

### Publication-Ready Outputs
- **High-quality figures** in multiple formats (PDF, PNG, SVG)
- **Customizable layouts** for different publication requirements
- **Integrated legends** and annotation systems
- **Scalable graphics** for presentations and manuscripts

## {octicon}`bookmark;0.85em` Citing CRESSENT

If you use CRESSENT in your work, please consider citing:

:::{card}
:link: [link to paper (still waiting)](https://www.biorxiv.org/)

**CRESSENT: A comprehensive toolkit for CRESS DNA virus analysis**
Pavan, R.R; Sullivan M.B.; Tisza, M. — *[Biorxiv]* ([2025]), DOI: [DOI].
:::

## {octicon}`people;0.85em` Related Tools

CRESSENT integrates with and builds upon several established bioinformatics tools:

- **MEME Suite** for de novo motif discovery
- **ViennaRNA** for RNA secondary structure prediction
- **IQ-TREE** for phylogenetic analysis
- **BLAST+** for sequence similarity searches
- **seqkit** for sequence manipulation and analysis

## {octicon}`question;0.85em` Ask a question or report a bug

If you want to ask a question about CRESSENT or report a problem, please create an issue in the [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/your-repo/cressent).

```{toctree}
:hidden:

self
```

```{toctree}
:caption: Using CRESSENT
:hidden:

installation
quickstart
faq
```

```{toctree}
:caption: Analysis Modules
:hidden:

align.md
phylogenetic_analysis
preprocessing
motif
secondary_structure
recombination
viz
```

