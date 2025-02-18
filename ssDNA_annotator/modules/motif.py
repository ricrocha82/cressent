#!/usr/bin/env python3
import argparse
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import logging

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
    """
    try:
        logging.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
        exit(1)

def find_pattern(input_fasta, seq_pattern, output_dir, table_name="pattern_positions.txt", remove_gaps=False):
    """
    Use seqkit to find the motif positions in the input FASTA file.
    The output table is saved in output_dir/table_name.
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
    with open(out2, "w") as seq2_file:
        seq2_file.write("\n".join(seq2_fasta))
    logging.info(f"Split sequences saved to: {out1} and {out2}")
    return out1, out2

def main():
    parser = argparse.ArgumentParser(
        description="Motif module to find a motif in sequences and optionally split them based on the motif position."
    )
    parser.add_argument("-i", "--input_fasta", required=True,
                        help="Input FASTA file with sequences (absolute or relative path).")
    parser.add_argument("-d", "--directory", default=".",
                        help="Directory for saving outputs and log files.")
    parser.add_argument("-p", "--pattern", required=True,
                        help="Sequence pattern (regex) for motif searching.")
    parser.add_argument("--split", action="store_true",
                        help="If set, the sequences will be split at the motif position.")
    parser.add_argument("-n", "--table_name", default="pattern_positions.txt",
                        help="Name of the file that will store motif positions (default: pattern_positions.txt).")
    parser.add_argument("--remove-gaps", action="store_true",
                        help="If set, removes gaps ('-') before searching for motifs.")
    
    args = parser.parse_args()
    output_dir = args.directory

    # Ensure the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory '{output_dir}' is ready.")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        exit(1)

    # Setup logging using the provided function
    setup_logging(output_dir)
    logging.info("Starting motif module.")
    
    # Determine the input FASTA full path (if a relative path is provided, assume it is in output_dir)
    input_fasta = args.input_fasta
    if not os.path.isabs(input_fasta):
        input_fasta = os.path.join(output_dir, input_fasta)
    logging.info(f"Using input FASTA: {input_fasta}")
    logging.info(f"Searching for motif with pattern: {args.pattern}")

    # Run motif finding using seqkit
    seq_table = find_pattern(input_fasta, args.pattern, output_dir, table_name=args.table_name, remove_gaps=args.remove_gaps)
    logging.info(f"Motif positions saved in: {seq_table}")

    # Optionally split sequences if --split flag is provided
    if args.split:
        prefix_out = os.path.join(output_dir, "split_sequences")
        split_fasta(seq_table, input_fasta, prefix_out)

    logging.info("Motif module finished.")

if __name__ == "__main__":
    main()


# python ./motif.py -i sequences.fasta -d /path/to/output -p "[GA].{4}GK[TS]" --split
