#!/usr/bin/env python3

import os
import sys
import importlib
import click
from pathlib import Path

# Get the absolute path to the directory containing this script
SCRIPT_DIR = Path(__file__).resolve().parent
MODULES_DIR = SCRIPT_DIR / "modules"

# Add the parent directory to the Python path to help import modules
sys.path.insert(0, str(SCRIPT_DIR.parent))

@click.group()
def cli():
    """
    Cressent: A comprehensive toolkit for ssDNA virus analysis
    
    Run a specific module by using: cressent [MODULE] [OPTIONS]
    
    To see help for a specific module, use: cressent [MODULE] --help
    """
    pass

# Build tree module
@cli.command(name="build_tree")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file with sequences.")
@click.option("-o", "--output", required=True, default = ".",help="Path to the output directory (Default: working directory)")
@click.option("-B", "--bootstrap", default=1000, type=int, help="Number of bootstrap iterations (default: 1000)")
@click.option("-t", "--threads", default="AUTO", help="Number of threads to use (default: AUTO)")
@click.option("-m", "--model", default="MFP", help="Substitution models (default: MFP - ModelFinder)")
@click.option("--keep_names", is_flag=True, help="Keep only the first word of sequence IDs, otherwise it replaces space with _")
@click.option("--extra_args", multiple=True, help="Extra arguments to pass directly to IQ-TREE (can be specified multiple times)")
def build_tree(input_fasta, output, bootstrap, threads, model, extra_args, keep_names):
    """Build phylogenetic tree using IQ-TREE."""
    try:
        from cressent.modules.build_tree import main as build_tree_main
        sys.argv = [sys.argv[0]] + ['--input_fasta', input_fasta, 
                                    '--output', output,
                                    '--bootstrap', str(bootstrap),
                                    '--threads', str(threads),
                                    '--model', model]
        # if extra_args:
        #     sys.argv.extend(['--extra_args', extra_args])
        if extra_args:
            # First add the flag
            sys.argv.append('--extra_args')
            # Then add each argument as a separate item
            sys.argv.extend(extra_args)
        
        if keep_names:
            sys.argv.append('--keep_names')

        build_tree_main()
    except Exception as e:
        click.echo(f"Error running build_tree module: {e}", err=True)
        sys.exit(1)

# Align module
@cli.command(name="align")
@click.option("-t", "--threads", type=int, default=1, help="Number of threads")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file with sequences")
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("--mafft_ep", type=float, default=0.123, help="Alignment length for MAFFT (default: 0.123)")
@click.option("--gap_threshold", type=float, default=0.2, help="Gap threshold for TrimAl (default: 0.2)")
@click.option("--db_family", help="List of family names for specific families, 'all' to use all the database, or 'custom' to use a custom AA file")
@click.option("--db_path", help="Path to the database FASTA files")
@click.option("--protein_type", type=click.Choice(['reps', 'caps']), help="Specify protein type (Rep or Cap) for database files")
@click.option("--custom_aa", help="Path to custom AA fasta file for alignment")
def align(threads, input_fasta, output, mafft_ep, gap_threshold, db_family, db_path, protein_type, custom_aa):
    """Pipeline for sequence alignment and trimming."""
    try:
        from cressent.modules.align import main as align_main
        sys.argv = [sys.argv[0]]
        if threads:
            sys.argv.extend(['--threads', str(threads)])
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if output:
            sys.argv.extend(['--output', output])
        if mafft_ep:
            sys.argv.extend(['--mafft_ep', str(mafft_ep)])
        if gap_threshold:
            sys.argv.extend(['--gap_threshold', str(gap_threshold)])
        if db_family:
            # Split by whitespace
            families = db_family.split()
            if families:
                sys.argv.append('--db_family')
                sys.argv.extend(families)
        if db_path:
            sys.argv.extend(['--db_path', db_path])
        if protein_type:
            sys.argv.extend(['--protein_type', protein_type])
        if custom_aa:
            sys.argv.extend(['--custom_aa', custom_aa])
        align_main()
    except Exception as e:
        click.echo(f"Error running align module: {e}", err=True)
        sys.exit(1)

# Cluster module
@cli.command(name="cluster")
@click.option("-i", "--input_fasta", required=True, help="Path to input FASTA file")
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("-t", "--threads", type=int, default=1, help="Number of threads for BLAST (Default = 1)")
@click.option("--min_ani", type=float, default=95.0, help="Minimum average identity for clustering (Default = 95.0)")
@click.option("--min_tcov", type=float, default=85.0, help="Minimum target coverage (Default = 85.0)")
@click.option("--min_qcov", type=float, default=0.0, help="Minimum query coverage (Default = 0.0)")
@click.option("--keep_names", is_flag=True, help="Keep only the first word of sequence IDs, otherwise replaces space with _")
def cluster(input_fasta, output, threads, min_ani, min_tcov, min_qcov, keep_names):
    """Sequence clustering using BLAST, anicalc, and aniclust."""
    try:
        from cressent.modules.cluster import main as cluster_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if output:
            sys.argv.extend(['--output', output])
        if threads:
            sys.argv.extend(['--threads', str(threads)])
        if min_ani:
            sys.argv.extend(['--min_ani', str(min_ani)])
        if min_tcov:
            sys.argv.extend(['--min_tcov', str(min_tcov)])
        if min_qcov:
            sys.argv.extend(['--min_qcov', str(min_qcov)])
        if keep_names:
            sys.argv.append('--keep_names')
        cluster_main()
    except Exception as e:
        click.echo(f"Error running cluster module: {e}", err=True)
        sys.exit(1)

# Motif discovery module
@cli.command(name="motif_disc")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file")
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory) (used for MEME and generated files)")
@click.option("-nmotifs", type=int, default=1, help="Number of motifs to find (Default = 1)")
@click.option("-minw", type=int, default=5, help="Minimum motif width (Default = 5)")
@click.option("-maxw", type=int, default=10, help="Maximum motif width (Default = 10)")
@click.option("--meme_extra", multiple=True, help="Additional MEME arguments (list format)")
@click.option("--scanprosite", is_flag=True, help='Run ScanProsite')
def motif_disc(input_fasta, output, nmotifs, minw, maxw, meme_extra, scanprosite):
    """Discover de novo motifs using MEME."""
    try:
        from cressent.modules.motif_disc import main as motif_disc_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if output:
            sys.argv.extend(['--output', output])
        if nmotifs:
            sys.argv.extend(['-nmotifs', str(nmotifs)])
        if minw:
            sys.argv.extend(['-minw', str(minw)])
        if maxw:
            sys.argv.extend(['-maxw', str(maxw)])
        if meme_extra:
            # First add the flag
            sys.argv.append('--meme_extra')
            # Then add each argument as a separate item
            sys.argv.extend(meme_extra)
        if scanprosite:
            sys.argv.append('--scanprosite')
        motif_disc_main()
    except Exception as e:
        click.echo(f"Error running motif_disc module: {e}", err=True)
        sys.exit(1)

# Motif gene map viz
@cli.command(name="motif_map_viz")
@click.option("-f","--file", help="Input file from scanprosite or motif analysis", required=True)
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("--format", type=click.Choice(['prosite', 'motif_table', 'auto']), default='auto', 
                help="Input format (default: auto-detect)")
def motif_map_viz(file, output, format):
    """Module to plot motif analysis results (ScanProsite or MEME motifs)"""
    try:
        from cressent.modules.motif_map_viz import main as motif_map_viz_main
        
        # Set up sys.argv to match argparse expectations
        sys.argv = [sys.argv[0]]
        sys.argv.extend(['-f', file])
        sys.argv.extend(['-o', output])
        sys.argv.extend(['--format', format])
        
        motif_map_viz_main()
        
    except Exception as e:
        click.echo(f"Error running motif_map_viz module: {e}", err=True)
        sys.exit(1)

# Motif module
@cli.command(name="motif")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file with sequences.")
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("-p", "--pattern", required=True, help="Sequence pattern (regex) for motif searching.")
@click.option("-n", "--table_name", default="pattern_positions.txt", help="Name of the file that will store motif positions (Default: pattern_positions.txt)")
@click.option("--remove-gaps", is_flag=True, help="If set, removes gaps ('-') before searching for motifs.")
@click.option("--split-sequences", is_flag=True, help="If set, the sequences will be split at the motif position.")
@click.option("--generate-logo", is_flag=True, help="If set, generate a sequence logo from the motif results.")
@click.option("--logo-name", default="sequence_logo.pdf", help="Name of the sequence logo PDF file (Default: sequence_logo.pdf)")
@click.option("--plot-title", default="sequence_logo", help="Title of the Sequence Logo (Default: sequence_logo)")
@click.option("--width", default=10, type=float, help="Width of the sequence logo PDF file (Default = 10)")
@click.option("--height", default=10, type=float, help="Height of the sequence logo PDF file (Default = 10)")
@click.option("--split-logo", is_flag=True, help="If set, the sequence logo will be split by group label.")
@click.option("--metadata", help="Path to metadata file containing group labels.")
@click.option("--ncol", type=int, help="Number of columns when splitting the sequence logo.")
@click.option("--group-label", help="Column name in metadata for grouping sequences.")
def motif(input_fasta, output, pattern, table_name, remove_gaps, split_sequences, 
          generate_logo, logo_name, plot_title, width, height, split_logo, 
          metadata, ncol, group_label):
    """Combined module for motif finding and sequence logo generation."""
    try:
        from cressent.modules.motif import main as motif_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if output:
            sys.argv.extend(['--output', output])
        if pattern:
            sys.argv.extend(['--pattern', pattern])
        if table_name:
            sys.argv.extend(['--table_name', table_name])
        if remove_gaps:
            sys.argv.append('--remove-gaps')
        if split_sequences:
            sys.argv.append('--split-sequences')
        if generate_logo:
            sys.argv.append('--generate-logo')
        if logo_name:
            sys.argv.extend(['--logo-name', logo_name])
        if plot_title:
            sys.argv.extend(['--plot-title', plot_title])
        if width:
            sys.argv.extend(['--width', str(width)])
        if height:
            sys.argv.extend(['--height', str(height)])
        if split_logo:
            sys.argv.append('--split-logo')
        if metadata:
            sys.argv.extend(['--metadata', metadata])
        if ncol:
            sys.argv.extend(['--ncol', str(ncol)])
        if group_label:
            sys.argv.extend(['--group-label', group_label])
        motif_main()
    except Exception as e:
        click.echo(f"Error running motif module: {e}", err=True)
        sys.exit(1)

# Gene map module
@cli.command(name="gene_map")
@click.option('-i', '--input', required=True, help='Input file path for the motif table CSV')
@click.option('-o', '--output', default='./', help="Path to the output directory (Default: working directory)")
@click.option('--filename', default='gene_motif.pdf', help='Output filename for the generated plot (default: gene_motif.pdf)')
@click.option('--height', type=float, default=10, help='Height of the output plot in inches (default: 10)')
@click.option('--width', type=float, default=10, help='Width of the output plot in inches (default: 10)')
@click.option('-t', '--title', help='Title for the plot (optional)')

def gene_map(input, output, filename, height, width, title):
    """Generate gene arrow plots from motif data using R."""
    try:
        from cressent.modules.gene_map import main as gene_map_main
        sys.argv = [sys.argv[0]]
        if input:
            sys.argv.extend(['--input', input])
        if output:
            sys.argv.extend(['--output', output])
        if filename:
            sys.argv.extend(['--filename', str(filename)])
        if height:
            sys.argv.extend(['--height', str(height)])
        if width:
            sys.argv.extend(['--width', str(width)])
        if title:
            sys.argv.extend(['--title', title])
        gene_map_main()
    except Exception as e:
        click.echo(f"Error running gene_map module: {e}", err=True)
        sys.exit(1)

# Plot tree module
@cli.command(name="plot_tree")
@click.option("-t", "--tree", help="Input tree file (Newick format)")
@click.option("--dist_matrix", help="Use distance matrix method for tree construction")
@click.option("-o", "--output", required=True, default=".",help="Path to the output directory (Default: working directory)")
@click.option("--metadata_1", help="Optional CSV metadata file")
@click.option("--metadata_2", help="Optional TSV name table file")
@click.option("--alignment", help="Optional alignment file (FASTA) to include in the plot")
@click.option("--layout", default="rectangular", help="Tree layout (e.g., rectangular, circular, unrooted) Default: rectangular")
@click.option("--branch_length", default="branch.length", help="Branch length parameter for ggtree (Default: branch.length)")
@click.option("--open_angle", default=0, type=float, help="Open angle for circular/unrooted layouts (default = 0)")
@click.option("--offset", default=0.14, type=float, help="Tip label offset (default = 0)")
@click.option("--tip_label", default="family", help="Column name to use as the tip label (default = family)")
@click.option("--color", default=True, help="Color tree by group (requires metadata) (default = True)")
@click.option("--fig_width", type=float, default=7, help="Figure width (ggsave) (default = 7)")
@click.option("--fig_height", type=float, default=7, help="Figure height (ggsave) (default = 7)")
@click.option("--plot_tips",  type=click.BOOL, default=True, help="Include tip labels in the plot (default = True)")
@click.option("--plot_name", default="tree_plot.pdf", help="Name of the output plot file (default: tree_plot.pdf)")
def plot_tree(tree, dist_matrix, output, metadata_1, metadata_2, alignment, layout, 
                branch_length, open_angle, offset, tip_label, color, fig_width, 
                fig_height, plot_tips, plot_name):
    """Plot phylogenetic trees using ggtree."""
    try:
        from cressent.modules.plot_tree import main as plot_tree_main
        sys.argv = [sys.argv[0]]
        if tree:
            sys.argv.extend(['--tree', tree])
        if dist_matrix:
            sys.argv.extend(['--dist_matrix', dist_matrix])
        if output:
            sys.argv.extend(['--output', output])
        if metadata_1:
            sys.argv.extend(['--metadata_1', metadata_1])
        if metadata_2:
            sys.argv.extend(['--metadata_2', metadata_2])
        if alignment:
            sys.argv.extend(['--alignment', alignment])
        if layout:
            sys.argv.extend(['--layout', layout])
        if branch_length:
            sys.argv.extend(['--branch_length', branch_length])
        if open_angle is not None:
            sys.argv.extend(['--open_angle', str(open_angle)])
        if offset is not None:
            sys.argv.extend(['--offset', str(offset)])
        if tip_label:
            sys.argv.extend(['--tip_label', tip_label])
        if color is False:
            sys.argv.extend(['--color', 'FALSE'])
        else:
            sys.argv.extend(['--color', 'TRUE'])
        if fig_width:
            sys.argv.extend(['--fig_width', str(fig_width)])
        if fig_height:
            sys.argv.extend(['--fig_height', str(fig_height)])
        if plot_tips is False:  # Explicitly check for False
            sys.argv.extend(['--plot_tips', 'FALSE'])
        else:
            sys.argv.extend(['--plot_tips', 'TRUE'])
        if plot_name:
            sys.argv.extend(['--plot_name', plot_name])
        plot_tree_main()
    except Exception as e:
        click.echo(f"Error running plot_tree module: {e}", err=True)
        sys.exit(1)

# Tanglegram module
@cli.command(name="tanglegram")
@click.option("--tree1", required=True, help="Path to the first tree file.")
@click.option("--tree2", required=True, help="Path to the second tree file.")
@click.option("--label1", required=True, help="Label for the first tree in the tanglegram.")
@click.option("--label2", required=True, help="Label for the second tree in the tanglegram.")
@click.option("-o", "--output", required=True, default=".",help="Path to the output directory (Default: working directory)")
@click.option("--name_tanglegram", default="tanglegram.pdf", help="Name of the tanglegram PDF file (default: tanglegram.pdf)")
@click.option("--width", default=20, type=float, help="Width of the tanglegram (default = 20)")
@click.option("--height", default=11, type=float, help="Height of the tanglegram (default = 11)")
@click.option("--lab_cex", default=1.5, type=float, help="cex size of the labels (default = 1.5)")
def tanglegram(tree1, tree2, label1, label2, output, name_tanglegram, width, height, lab_cex):
    """Generate a tanglegram from two phylogenetic trees."""
    try:
        from cressent.modules.tanglegram import main as tanglegram_main
        sys.argv = [sys.argv[0]]
        if tree1:
            sys.argv.extend(['--tree1', tree1])
        if tree2:
            sys.argv.extend(['--tree2', tree2])
        if label1:
            sys.argv.extend(['--label1', label1])
        if label2:
            sys.argv.extend(['--label2', label2])
        if output:
            sys.argv.extend(['--output', output])
        if name_tanglegram:
            sys.argv.extend(['--name_tanglegram', name_tanglegram])
        if width:
            sys.argv.extend(['--width', str(width)])
        if height:
            sys.argv.extend(['--height', str(height)])
        if lab_cex:
            sys.argv.extend(['--lab_cex', str(lab_cex)])
        tanglegram_main()
    except Exception as e:
        click.echo(f"Error running tanglegram module: {e}", err=True)
        sys.exit(1)

# Sequence Logo module
@cli.command(name="seq_logo")
@click.option("-i", "--input_fasta", help="Path to the fasta file.")
@click.option("-tb", "--seq_df", help="Path to the table produced by seqkit.")
@click.option("-o", "--output", default="./", help="Path to the output directory (Default: working directory)")
@click.option("-n", "--output_name", default="sequence_logo.pdf", help="Name of the sequence logo (default: sequence_logo.pdf)")
@click.option("--plot_title", default="sequence_logo", help="Title of the Sequence Logo (default: sequence_logo)")
@click.option("--width", default=10, type=float, help="Width of the sequence logo (default = 10)")
@click.option("--height", default=10, type=float, help="Height of the sequence logo (default = 10)")
@click.option("--split", is_flag=True, help="If set, the sequence logo will be split by group label (default = True)")
@click.option("--metadata", help="Path to metadata file containing group labels.")
@click.option("--ncol", type=int, help="Number of columns when splitting the sequence logo.")
@click.option("--group_label", help="Column name in metadata for grouping sequences.")
@click.option("--positions_per_row", default=50, type=int, help="Number of positions per row when creating multi-row plots (default: 50)")
@click.option("--max_positions_single_row", default=100, type=int, help="Maximum number of positions before automatically splitting into multiple rows (default: 100)")
@click.option("--method", type=click.Choice(["bits", "prob"]), default="prob", help="Method for ggseqlogo: 'bits' for information content or 'prob' for probability (default: prob)")
def seq_logo(input_fasta, seq_df, output, output_name, plot_title, width, height, split, metadata, ncol, group_label, positions_per_row, max_positions_single_row, method):
    """Generate a sequence logo from a FASTA file or sequence table."""
    try:
        from cressent.modules.seq_logo import main as seq_logo_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['-i', input_fasta])
        if seq_df:
            sys.argv.extend(['-tb', seq_df])
        if output:
            sys.argv.extend(['-o', output])
        if output_name:
            sys.argv.extend(['-n', output_name])
        if plot_title:
            sys.argv.extend(['--plot_title', plot_title])
        if width:
            sys.argv.extend(['--width', str(width)])
        if height:
            sys.argv.extend(['--height', str(height)])
        if split:
            sys.argv.append('--split')
        if metadata:
            sys.argv.extend(['--metadata', metadata])
        if ncol:
            sys.argv.extend(['--ncol', str(ncol)])
        if group_label:
            sys.argv.extend(['--group_label', group_label])
        if positions_per_row:
            sys.argv.extend(['--positions_per_row', str(positions_per_row)])
        if max_positions_single_row:
            sys.argv.extend(['--max_positions_single_row', str(max_positions_single_row)])
        if method:
            sys.argv.extend(['--method', method])
        seq_logo_main()
    except Exception as e:
        click.echo(f"Error running seq_logo module: {e}", err=True)
        sys.exit(1)

# Stem-loop finder module
@cli.command(name="sl_finder")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file")
@click.option("--gff_in", required=True, help="Input GFF/GTF file")
@click.option("--out_gff", required=True, help="Output GFF filename")
@click.option("-o","--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("--csv_out", help="Output CSV filename")
@click.option("--motif", default="nantantan", help="Conserved motif (default = nantantan)")
@click.option("--family", type=click.Choice(['geminiviridae', 'genomoviridae', 'smacoviridae', 
                                            'cycloviridae', 'circoviridae', 'general']), 
             help="CRESS viral family")
@click.option("--idealstemlen", "-s", type=int, default=11, help="Ideal stem length (default = 11)")
@click.option("--ideallooplen", "-l", type=int, default=11, help="Ideal loop length (default = 11)")
@click.option("--frame", "-f", type=int, default=15, help="Bases around motif for folding (default = 15)")
def sl_finder(input_fasta, gff_in, out_gff, output, csv_out, motif, family, 
              idealstemlen, ideallooplen, frame):
    """A module for putative stem-loop annotation."""
    try:
        from cressent.modules.sl_finder import main as sl_finder_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['-i', input_fasta])
        if gff_in:
            sys.argv.extend(['--gff_in', gff_in])
        if out_gff:
            sys.argv.extend(['--out_gff', out_gff])
        if output:
            sys.argv.extend(['--output', output])
        if csv_out:
            sys.argv.extend(['--csv_out', csv_out])
        if motif:
            sys.argv.extend(['--motif', motif])
        if family:
            sys.argv.extend(['--family', family])
        if idealstemlen:
            sys.argv.extend(['--idealstemlen', str(idealstemlen)])
        if ideallooplen:
            sys.argv.extend(['--ideallooplen', str(ideallooplen)])
        if frame:
            sys.argv.extend(['--frame', str(frame)])
        # Call the main function, but we need to create it first in the sl_finder.py file
        sl_finder_main()
    except Exception as e:
        click.echo(f"Error running sl_finder module: {e}", err=True)
        sys.exit(1)

# GC Patchiness module (heatmap)
@cli.command(name="gc_ht")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file containing sequences.")
@click.option("-o","--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("--window_size", type=int, default=30, help="Sliding window size for GC content calculation (default: 30).")
@click.option("--step_size", type=int, default=5, help="Step size for GC calculation (default: 5).")
@click.option("--xticklabels", type=int, default=10, help="Interval for x-axis tick labels (default: None).")
@click.option("--fig_width", type=int, default=12, help="Figure width in inches (default: 12).")
@click.option("--fig_height", type=int, default=6, help="Figure height in inches (default: 6).")
@click.option("--output_name", default="gc_heatmap.pdf", help="Name of output image file with extension (default: gc_heatmap.png).")
def gc_ht(input_fasta, output, window_size, step_size, xticklabels, fig_width, fig_height, output_name):
    """Generate a GC content heatmap from a FASTA file."""
    try:
        from cressent.modules.gc_ht import main as gc_ht_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if output:
            sys.argv.extend(['--output', output])
        if window_size:
            sys.argv.extend(['--window_size', str(window_size)])
        if step_size:
            sys.argv.extend(['--step_size', str(step_size)])
        if xticklabels:
            sys.argv.extend(['--xticklabels', str(xticklabels)])
        if fig_width:
            sys.argv.extend(['--fig_width', str(fig_width)])
        if fig_height:
            sys.argv.extend(['--fig_height', str(fig_height)])
        if output_name:
            sys.argv.extend(['--output_name', output_name])
        # Call the main function
        gc_ht_main()
    except Exception as e:
        click.echo(f"Error running gc_ht module: {e}", err=True)
        sys.exit(1)

# Recombination module
@cli.command(name="recombination")
@click.option('-i', '--input_fasta', required=True, help='Input alignment file in FASTA format')
@click.option('-o', '--output', default='.', help='Path to the output directory (Default: working directory)')
@click.option('-f', '--output_file', required=True, help='Output file for results (CSV format)')
@click.option('-c', '--config', help='Configuration file in INI format for OpenRDP parameters')
@click.option('-rdp', is_flag=True, help='Run RDP method')
@click.option('-threeseq', is_flag=True, help='Run 3Seq method')
@click.option('-geneconv', is_flag=True, help='Run GENECONV method')
@click.option('-maxchi', is_flag=True, help='Run MaxChi method')
@click.option('-chimaera', is_flag=True, help='Run Chimaera method')
@click.option('-bootscan', is_flag=True, help='Run Bootscan method')
@click.option('-siscan', is_flag=True, help='Run Siscan method')
@click.option('-all', is_flag=True, help='Run all methods')
@click.option('-quiet', is_flag=True, help='Suppress console output')
@click.option('-verbose', is_flag=True, help='Enable verbose logging')
def recombination(input_fasta, output, output_file, config, rdp, threeseq, geneconv, maxchi, 
                  chimaera, bootscan, siscan, all, quiet, verbose):
    """Detect recombination events in ssDNA virus sequences."""
    try:
        from cressent.modules.recombination import main as recombination_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['-i', input_fasta])
        if output:
            sys.argv.extend(['-o', output])
        if output_file:
            sys.argv.extend(['-f', output_file])
        if config:
            sys.argv.extend(['-c', config])
        if rdp:
            sys.argv.append('-rdp')
        if threeseq:
            sys.argv.append('-threeseq')
        if geneconv:
            sys.argv.append('-geneconv')
        if maxchi:
            sys.argv.append('-maxchi')
        if chimaera:
            sys.argv.append('-chimaera')
        if bootscan:
            sys.argv.append('-bootscan')
        if siscan:
            sys.argv.append('-siscan')
        if all:
            sys.argv.append('-all')
        if quiet:
            sys.argv.append('-quiet')
        if verbose:
            sys.argv.append('-verbose')
        recombination_main()
    except Exception as e:
        click.echo(f"Error running recombination module: {e}", err=True)
        sys.exit(1)

# Run CRUISE module
@cli.command(name="run_cruise")
@click.option("-o","--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("-i","--input_fasta", required=True, help="Path to input FASTA file with all sequences")
@click.option("--inputGFF", required=True, help="Path to associated input GFF file")
@click.option("--outputGFF", required=True, default="finaloutput.gff",help="Path for output GFF file (default: finaloutput.gff)")
@click.option("--outputAnnotations", default="2 CRUISE", help="Identifiers to selectively preserve annotations (default: 2 CRUISE)")
@click.option("--minLength", type=int, default=5, help="Minimum iteron length (default = 5)")
@click.option("--maxLength", type=int, default=12, help="Maximum iteron length (default = 12)")
@click.option("--range", type=int, default=65, help="Number of base pairs around nona to search (default = 65)")
@click.option("--rank", type=bool, default=True, help="Use ranking system (default = True)")
@click.option("--numberTopIterons", type=int, default=5, help="The number of iterons returned in rank order (default = 5)")
@click.option("--maxScore", type=int, default=40, help="Maximum score allowed for iterons if rank = False (default = 40)")
@click.option("--wiggle", type=int, default=5, help="Max difference between iteron length and distance (default = 5)")
@click.option("--goodLength", type=int, default=11, help="The highest favorable iteron length (default = 11)")
@click.option("--doStemLoop", type=bool, default=True, help="Whether to annotate stem-loop repeats (default = True)")
@click.option("--doKnownIterons", type=bool, default=True, help="Whether to annotate known iterons (default = True)")
@click.option("--maxDist", type=int, default=20, help="Maximum allowed distance between iterons (default = 20)")
@click.option("--bestDist", type=int, default=10, help="Optimal maximum distance between iterons (default = 10)")
@click.option("--scoreRange", type=int, default=50, help="Score range between outputted candidates (default = 50)")
def run_cruise(output, input_fasta, inputgff, outputgff, outputannotations,
              minlength, maxlength, range, rank, numbertopiterons, maxscore,
              wiggle, goodlength, dostemloop, doknowniterons, maxdist, bestdist, scorerange):
    """Search for iterons around CRESS stem-loops in GFF files."""
    try:
        from cressent.modules.run_cruise import main as run_cruise_main
        sys.argv = [sys.argv[0]]
        if output:
            sys.argv.extend(['--output', outputdir])
        if input_fasta:
            sys.argv.extend(['--input_fasta', input_fasta])
        if inputgff:
            sys.argv.extend(['--inputGFF', inputgff])
        if outputgff:
            sys.argv.extend(['--outputGFF', outputgff])
        if outputannotations:
            sys.argv.extend(['--outputAnnotations', outputannotations])
        if minlength:
            sys.argv.extend(['--minLength', str(minlength)])
        if maxlength:
            sys.argv.extend(['--maxLength', str(maxlength)])
        if range:
            sys.argv.extend(['--range', str(range)])
        if rank is not None:
            sys.argv.extend(['--rank', str(rank).lower()])
        if numbertopiterons:
            sys.argv.extend(['--numberTopIterons', str(numbertopiterons)])
        if maxscore:
            sys.argv.extend(['--maxScore', str(maxscore)])
        if wiggle:
            sys.argv.extend(['--wiggle', str(wiggle)])
        if goodlength:
            sys.argv.extend(['--goodLength', str(goodlength)])
        if dostemloop is not None:
            sys.argv.extend(['--doStemLoop', str(dostemloop).lower()])
        if doknowniterons is not None:
            sys.argv.extend(['--doKnownIterons', str(doknowniterons).lower()])
        if maxdist:
            sys.argv.extend(['--maxDist', str(maxdist)])
        if bestdist:
            sys.argv.extend(['--bestDist', str(bestdist)])
        if scorerange:
            sys.argv.extend(['--scoreRange', str(scorerange)])
        run_cruise_main()
    except Exception as e:
        click.echo(f"Error running run_cruise module: {e}", err=True)
        sys.exit(1)

# Detect contamination module
@cli.command(name="detect_contamination")
@click.option("-i", "--input_fasta", required=True, help="Input FASTA file")
@click.option("--db", required=True, help="Contaminant database FASTA file")
@click.option("-o", "--output", required=True, help="Path to the output directory (Default: working directory)")
@click.option("--output-name", default="clean_sequences", help="Base name for output files (default: clean_sequences)")
@click.option("--seq-type", type=click.Choice(['nucl', 'prot']), help="Sequence type (auto-detect if not specified)")
@click.option("--evalue", type=float, default=1e-10, help="BLAST E-value threshold (default: 1e-10)")
@click.option("--identity", type=float, default=90.0, help="Minimum percent identity to consider a match (default: 90.0)")
@click.option("--coverage", type=float, default=50.0, help="Minimum query coverage to consider a match (default: 50.0)")
@click.option("-t","--threads", type=int, default=1, help="Number of CPU threads for BLAST (default: 1)")
@click.option("--keep-temp", is_flag=True, help="Keep temporary BLAST output files")
def detect_contamination(input_fasta, db, output, output_name, seq_type, evalue, identity, 
                        coverage, threads, keep_temp):
    """Filter viral contaminants from sequence data (nucleotide or protein)."""
    try:
        from cressent.modules.detect_contamination import main as detect_contamination_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['-i', input_fasta])
        if db:
            sys.argv.extend(['--db', db])
        if output:
            sys.argv.extend(['-o', output])
        if output_name:
            sys.argv.extend(['--output-name', output_name])
        if seq_type:
            sys.argv.extend(['--seq-type', seq_type])
        if evalue:
            sys.argv.extend(['--evalue', str(evalue)])
        if identity:
            sys.argv.extend(['--identity', str(identity)])
        if coverage:
            sys.argv.extend(['--coverage', str(coverage)])
        if threads:
            sys.argv.extend(['-t', str(threads)])
        if keep_temp:
            sys.argv.append('--keep-temp')
        detect_contamination_main()
    except Exception as e:
        click.echo(f"Error running detect_contamination module: {e}", err=True)
        sys.exit(1)

# Adjust sequence start position module
@cli.command(name="adjust_seq")
@click.option("-i", "--input_fasta", required=True, help="Path to the input FASTA file.")
@click.option("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
@click.option("-m", "--motif", default="TAGTATTAC", help="Motif to adjust sequences to start with (default: TAGTATTAC).")
def adjust_seq(input_fasta, output, motif):
    """Adjust sequences in a FASTA file to start with a specified motif."""
    try:
        from cressent.modules.adjust_seq import main as adjust_seq_main
        sys.argv = [sys.argv[0]]
        if input_fasta:
            sys.argv.extend(['-i', input_fasta])
        if output:
            sys.argv.extend(['-o', output])
        if motif:
            sys.argv.extend(['-m', motif])
        adjust_seq_main()
    except Exception as e:
        click.echo(f"Error running adjust_seq module: {e}", err=True)
        sys.exit(1)

# Build contaminant database module
@cli.command(name="build_contaminant_db")
@click.option("--accession-csv", required=True, help="CSV file with accessions (must have 'accession' column)")
@click.option("-o","--output", required=True, help="Path to the output directory (Default: working directory)")
@click.option("--output-name", default="contaminant_db", help="Base name for output files (default: contaminant_db)")
@click.option("--email", default="user@example.com", help="Email for NCBI Entrez queries (not required, default: user@example.com)")
@click.option("--batch-size", type=int, default=10, help="Maximum number of sequences to download in each batch (default = 10)")
def build_contaminant_db(accession_csv, output, output_name, email, batch_size):
    """Build a viral contaminant database for decontamination pipelines."""
    try:
        from cressent.modules.build_contaminant_db import main as build_contaminant_db_main
        sys.argv = [sys.argv[0]]
        if accession_csv:
            sys.argv.extend(['--accession-csv', accession_csv])
        if output:
            sys.argv.extend(['--output', output])
        if output_name:
            sys.argv.extend(['--output-name', output_name])
        if email:
            sys.argv.extend(['--email', email])
        if batch_size:
            sys.argv.extend(['--batch-size', str(batch_size)])
        build_contaminant_db_main()  # Actually call the main function
    except Exception as e:
        click.echo(f"Error running build_contaminant_db module: {e}", err=True)
        sys.exit(1)


# Build my own DB
@cli.command(name="db_builder")
@click.option('-t', '--taxonomy-file', required=True, help='Path to taxonomy_accession_number.csv file')
@click.option('-l', '--taxonomy-level', required=True,
              type=click.Choice(['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 
                               'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 
                               'Genus', 'Subgenus', 'Species']),
              help='Taxonomy level to use for selection')
@click.option('-s', '--selected-taxonomies', help='Selected taxonomies (if not provided, will list available options)')
@click.option('-o', '--output-dir', default=".", help='Output directory for database')
@click.option('-e', '--email', default="my_email@gmail.com", help='Email for NCBI Entrez')
@click.option('--threads', type=int, default=8, help='Number of threads to use (default: 8)')
@click.option('--cd-hit-identity', type=float, default=0.95, help='CD-HIT identity threshold (default: 0.95)')
@click.option('--mcl-inflation', type=float, default=1.5, help='MCL inflation parameter (default: 1.5)')
def db_builder(taxonomy_file, taxonomy_level, selected_taxonomies, output_dir, email, threads, cd_hit_identity, mcl_inflation):
    """ Build taxonomy-based database for ssDNA tool """
    try:
        from cressent.modules.db_builder import main as db_builder_main
        
        # Construct sys.argv for the db_builder module
        sys.argv = ['db_builder.py']
        sys.argv.extend(['--taxonomy-file', taxonomy_file])
        sys.argv.extend(['--taxonomy-level', taxonomy_level])
        
        if selected_taxonomies:
            # Split by whitespace and add to sys.argv
            taxa = selected_taxonomies.split()
            if taxa:
                sys.argv.extend(['--selected-taxonomies'] + taxa)
        
        sys.argv.extend(['--output-dir', output_dir])
        sys.argv.extend(['--email', email])
        sys.argv.extend(['--threads', str(threads)])
        sys.argv.extend(['--cd-hit-identity', str(cd_hit_identity)])
        sys.argv.extend(['--mcl-inflation', str(mcl_inflation)])
        
        db_builder_main()  # Actually call the main function
    except Exception as e:
        click.echo(f"Error running db_builder module: {e}", err=True)
        sys.exit(1)

# Dynamically load additional modules
def load_dynamic_modules():
    """
    Dynamically discover and load modules from the modules directory.
    This allows new modules to be added without modifying this file.
    """
    # Skip modules that we've already added explicitly
    existing_commands = [cmd.name for cmd in cli.commands.values()]
    
    try:
        for file_path in MODULES_DIR.glob("*.py"):
            # Skip __init__.py and already loaded modules
            module_name = file_path.stem
            if module_name.startswith("__") or module_name in existing_commands:
                continue
                
            # Try to import the module to see if it has a main function
            try:
                spec = importlib.util.spec_from_file_location(module_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                # Check if the module has a main function
                if hasattr(module, "main"):
                    # Create a command for this module
                    click.echo(f"Dynamically loaded module: {module_name}")
                    # Note: This would be much more complex to implement properly
                    # as you'd need to dynamically create Click commands with proper arguments
            except Exception as e:
                # Skip modules that can't be loaded
                pass
    except Exception as e:
        # Fail silently if dynamic loading doesn't work
        pass

if __name__ == "__main__":
    # Load dynamic modules - commented out for now as it needs more work
    # load_dynamic_modules()
    cli()