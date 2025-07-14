#!/usr/bin/env python3
import argparse
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import logging
import sys
from pathlib import Path

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "motif.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def run_command(command, error_message):
    """
    Run a shell command and log its progress. Exits on error.
    
    Args:
        command: The command to run (either string or list)
        error_message: Message to log if the command fails
    """
    try:
        logging.info(f"Running command: {command}")
        if isinstance(command, list):
            subprocess.run(command, check=True)
        else:
            subprocess.run(command, shell=True, check=True)
        logging.info(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
        exit(1)

def validate_fasta(filename, output_dir=None):
    """
    Validate that a file is in FASTA format.
    
    Args:
        filename: Path to the FASTA file
        output_dir: Optional output directory to prepend to the filename
        
    Returns:
        The full path to the validated FASTA file
    """
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                logging.info("FASTA format validated.")

                return filename
            else:
                logging.error(f"Error: Input file {filename} is not in the FASTA format.")
                sys.exit("Error: Input file is not in the FASTA format.\n")
    except Exception as e:
        logging.error(f"Error validating FASTA file: {e}")
        sys.exit(f"Error: Could not open or read file {filename}.\n")

def find_pattern(input_fasta, seq_pattern, output_dir, table_name="pattern_positions.txt", remove_gaps=False):
    """
    Use seqkit to find the motif positions in the input FASTA file.
    The output table is saved in output_dir/table_name.
    
    Args:
        input_fasta: Path to the input FASTA file
        seq_pattern: Sequence pattern (regex) to search for
        output_dir: Directory where the output table will be saved
        table_name: Name of the output table file
        remove_gaps: Whether to remove gaps from sequences before searching
        
    Returns:
        Path to the created table file
    """
    seq_table = os.path.join(output_dir, table_name)

    if remove_gaps:
        ungapped_fasta = os.path.join(output_dir, "ungapped_sequences.fasta")
        run_command(f"seqkit seq --remove-gaps {input_fasta} > {ungapped_fasta}", "Error removing gaps")
        input_fasta = ungapped_fasta  # Use the ungapped file

    cmd = f"seqkit locate --id-ncbi -i -r -p '{seq_pattern}' {input_fasta} > {seq_table}"
    run_command(cmd, "Error locating patterns with SeqKit")

    try:
        with open(seq_table, "r") as f:
            lines = f.readlines()

        # Check if file has more than just a header
        if len(lines) > 1:
            num_patterns = len(lines) - 1  # Subtract header line
        else:
            num_patterns = 0  # No motifs found

        logging.info(f"Total number of patterns found: {num_patterns}")
        
    except Exception as e:
        logging.error(f"Error counting patterns: {e}")
        return 0

    return seq_table

def split_fasta(table_file, fasta_file, prefix_out):
    """
    Split each sequence in fasta_file at the first motif occurrence.
    Two FASTA files are written with prefixes <prefix_out>_1.fasta and <prefix_out>_2.fasta.
    
    Args:
        table_file: Path to the table file with motif positions
        fasta_file: Path to the FASTA file with sequences
        prefix_out: Prefix for output files
        
    Returns:
        Tuple of (path to first FASTA file, path to second FASTA file)
    """

    logging.info(f"Reading motif table: {table_file}")
    split_table = pd.read_csv(table_file, sep="\t")
    logging.info(f"Reading FASTA file: {fasta_file}")
    fasta_sequences = {record.description: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    seq1_fasta = []
    seq2_fasta = []

    for _, row in split_table.iterrows():
        seq_id = row['seqID']
        start = int(row['start'])
        if seq_id in fasta_sequences:
            sequence = fasta_sequences[seq_id]
            # Split at the motif position:
            seq1_part = sequence[:start - 1]  # from beginning up to (start-1)
            seq2_part = sequence[start - 1:]  # from motif start onward
            seq1_fasta.append(f">{seq_id}\n{seq1_part}")
            seq2_fasta.append(f">{seq_id}\n{seq2_part}")
            logging.info(f"Sequence '{seq_id}' split at position {start}.")
        else:
            logging.warning(f"Sequence ID '{seq_id}' not found in the FASTA file.")

    out1 = f"{prefix_out}_1.fasta"
    out2 = f"{prefix_out}_2.fasta"
    with open(out1, "w") as seq1_file:
        seq1_file.write("\n".join(seq1_fasta))
        seq1_file.write("\n")
    with open(out2, "w") as seq2_file:
        seq2_file.write("\n".join(seq2_fasta))
        seq2_file.write("\n")
    logging.info(f"Split sequences saved to: {out1} and {out2}")
    return out1, out2

def generate_seqlogo(seq_table, output_dir, output_name="sequence_logo.pdf",
                     plot_title="sequence_logo", width=10, height=10, split=False,
                     metadata=None, ncol=None, group_label=None):
    """
    Generate a sequence logo using an R script.

    Args:
        fasta_file: Path to the FASTA file (or None if only using table)
        seq_table: Path to the table that seqkit generated
        output_dir: Directory to save the output
        output_name: Name of the output PDF file
        plot_title: Title of the Sequence Logo
        width: Width of the PDF file
        height: Height of the PDF file
        split: Whether to split by group labels
        metadata: Path to the metadata file for grouping
        ncol: Number of columns in the split plot
        group_label: Column name in metadata for grouping sequences
        
    Returns:
        Path to the generated logo file
    """
    # Path to the R script
    modules_dir = Path(__file__).resolve().parent
    r_script_path = modules_dir / "seq_logo.R"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Output directory set to: {output_dir}")

    # Build the command
    cmd = [
        "Rscript",
        str(r_script_path),
        f"seq_df={seq_table if seq_table else ''}",
        f"plot_title={plot_title}",
        f"output_dir={output_dir}",
        f"output_name={output_name}",
        f"height={height}",
        f"width={width}",
        f"split={'TRUE' if split else 'FALSE'}",
        f"metadata={metadata if metadata else ''}",
        f"ncol={ncol if ncol else ''}",
        f"group_label={group_label if group_label else ''}"
    ]

    # Run the R script
    run_command(cmd, "Error generating sequence logo")
    
    return os.path.join(output_dir, output_name)

def main():
    parser = argparse.ArgumentParser(
        description="Combined module for motif finding and sequence logo generation."
    )
    
    # Input file options
    parser.add_argument("-i", "--input_fasta", required=True,
                        help="Input FASTA file with sequences.")
    
    # Output options
    parser.add_argument("-o", "--output", default=".",
                        help="Directory for saving outputs and log files.")
    
    # Motif finding options
    parser.add_argument("-p", "--pattern", required=True,
                        help="Sequence pattern (regex) for motif searching.")
    parser.add_argument("-n", "--table_name", default="pattern_positions.txt",
                        help="Name of the file that will store motif positions (Default: pattern_positions.txt)")
    parser.add_argument("--remove-gaps", action="store_true",
                        help="If set, removes gaps ('-') before searching for motifs.")
    parser.add_argument("--split-sequences", action="store_true",
                        help="If set, the sequences will be split at the motif position.")

    # Logo generation options
    parser.add_argument("--generate-logo", action="store_true",
                        help="If set, generate a sequence logo from the motif results.")
    parser.add_argument("--logo-name", default="sequence_logo.pdf",
                        help="Name of the sequence logo PDF file (Default: sequence_logo.pdf)")
    parser.add_argument("--plot-title", default="sequence_logo",
                        help="Title of the Sequence Logo (Default: sequence_logo)")
    parser.add_argument("--width", default=10, type=float,
                        help="Width of the sequence logo PDF file (Default = 10)")
    parser.add_argument("--height", default=10, type=float,
                        help="Height of the sequence logo PDF file (Default = 10)")
    parser.add_argument("--split-logo", action="store_true",
                        help="If set, the sequence logo will be split by group label.")
    parser.add_argument("--metadata", 
                        help="Path to metadata file containing group labels.")
    parser.add_argument("--ncol", type=int,
                        help="Number of columns when splitting the sequence logo.")
    parser.add_argument("--group-label", 
                        help="Column name in metadata for grouping sequences.")
    
    args = parser.parse_args()

    # Ensure the output directory exists
    output_dir = args.output
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory '{output_dir}' is ready.")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        exit(1)

    # Setup logging using the provided function
    setup_logging(output_dir)
    logging.info("Starting motif module.")

    # Validate and process FASTA file
    input_fasta = validate_fasta(args.input_fasta, output_dir)
    logging.info(f"Using input FASTA: {input_fasta}")
    logging.info(f"Searching for motif with pattern: {args.pattern}")
    
    # if not os.path.isabs(input_fasta):
    #     input_fasta = os.path.join(output_dir, input_fasta)

    # Run motif finding using seqkit
    seq_table = find_pattern(input_fasta, args.pattern, output_dir, table_name=args.table_name, remove_gaps=args.remove_gaps)
    if not seq_table:
        logging.error("Motif finding failed or no motifs found.")
        exit(1)
    logging.info(f"Motif positions saved in: {seq_table}")

    # Optionally split sequences if --split-sequences flag is provided
    if args.split_sequences:
        prefix_out = os.path.join(output_dir, "split_sequences")
        split_files = split_fasta(seq_table, input_fasta, prefix_out)
        logging.info(f"Sequences split into: {split_files[0]} and {split_files[1]}")

    # Optionally generate sequence logo if --generate-logo flag is provided
    if args.generate_logo:
        # Validate inputs for logo generation with splitting
        if args.split_logo and (args.metadata is None or args.group_label is None or args.ncol is None):
            logging.error("When using --split-logo, you must provide --metadata, --group-label, and --ncol.")
            print("Error: When using --split-logo, you must provide --metadata, --group-label, and --ncol.")
            exit(1)
        
        logo_file = generate_seqlogo(
            seq_table=seq_table,
            output_dir=output_dir,
            output_name=args.logo_name,
            plot_title=args.plot_title,
            width=args.width,
            height=args.height,
            split=args.split_logo,
            metadata=args.metadata,
            ncol=args.ncol,
            group_label=args.group_label
        )
        logging.info(f"Sequence logo generated: {logo_file}")

    logging.info("Motif logo module finished.")

if __name__ == "__main__":
    main()


