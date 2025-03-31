import click
from ssDNA_annotator.modules.align import detect_splicing
from ssDNA_annotator.phylogeny.build_tree import build_phylogenetic_tree
from ssDNA_annotator.motifs.find_motifs import find_motifs
from ssDNA_annotator.decontamination.decontaminate import decontaminate
from ssDNA_annotator.visualization.make_figures import create_figures

@click.group()
def cli():
    """CLI for ssDNA Virus Annotation Tool"""
    pass

@cli.command()
@click.argument("input_file")
def splicing(input_file):
    """Detect splicing (introns and exons)"""
    detect_splicing(input_file)

@cli.command()
@click.argument("input_file")
def phylogeny(input_file):
    """Build phylogenetic tree using IQ-TREE"""
    build_phylogenetic_tree(input_file)

@cli.command()
@click.argument("input_file")
def motifs(input_file):
    """Find motifs using Prodigal and Prodigal-GV"""
    find_motifs(input_file)

@cli.command()
@click.argument("input_file")
@click.option("--database", required=True, help="Path to contamination database")
def decontamination(input_file, database):
    """Decontaminate sequences based on a database"""
    decontaminate(input_file, database)

@cli.command()
@click.argument("input_file")
def visualization(input_file):
    """Create genome visualization figures"""
    create_figures(input_file)

if __name__ == "__main__":
    cli()