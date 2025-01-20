#!/usr/bin/env python3

import argparse
import subprocess
import os
import re
import pandas as pd
from Bio import SeqIO


def run_command(command, error_message):
    """
    Helper function to run a shell command and handle errors.
    """
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        print(f"{error_message}: {e}")
        exit(1)

def find_and_split_by_pattern(input_fasta, seq_pattern, output_dir):
    """
    Find pattern positions and split sequences based on the pattern.
    """
    seq_table = os.path.join(output_dir, "pattern_positions.txt")
    prefix_out = os.path.join(output_dir, "split_sequences")

    # Find pattern positions using seqkit
    cmd = f"seqkit locate --id-ncbi -i -r -p '{seq_pattern}' {input_fasta} > {seq_table}"
    run_command(cmd, "Error locating patterns with SeqKit")

    # Split sequences into two domains
    split_fasta(seq_table, input_fasta, prefix_out)
    return f"{prefix_out}_1.fasta", f"{prefix_out}_2.fasta"


def split_fasta(table_file, fasta_file,prefix_out):
    """
    Split sequences in the FASTA file at specific pattern positions.
    """

    # Step 1: Load the table
    split_table = pd.read_csv(table_file, sep="\t")

    # Step 2: Read the FASTA file
    fasta_sequences = {record.description: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    # Step 3: Initialize outputs
    seq1_fasta = []
    seq2_fasta = []

    # Process each row in the table
    for _, row in split_table.iterrows():
        seq_id = row['seqID']
        start = int(row['start'])

        # Ensure the sequence exists in the FASTA file
        if seq_id in fasta_sequences:
            sequence = fasta_sequences[seq_id]
            
            # Split the sequence
            seq1_part = sequence[:start - 1]  # Up to (start - 1)
            seq2_part = sequence[start - 1:]      # From start onward

            # Append to respective outputs
            seq1_fasta.append(f">{seq_id}\n{seq1_part}")
            seq2_fasta.append(f">{seq_id}\n{seq2_part}")
        else:
            print(f"Warning: Sequence ID '{seq_id}' not found in FASTA file.")

    # Step 5: Write outputs to separate files
    with open(f"{prefix_out}_1.fasta", "w") as seq1_file:
        seq1_file.write("\n".join(seq1_fasta))

    with open(f"{prefix_out}_2.fasta", "w") as seq2_file:
        seq2_file.write("\n".join(seq2_fasta))

def main():
    parser = argparse.ArgumentParser(description="Pipeline to split a sequence based on determined patter (e.g., motif)")
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing input FASTA file and for saving outputs.")
    parser.add_argument("-p", "--pattern", required=True, help="Sequence pattern (regex) for splitting sequences.")

    args = parser.parse_args()
    input_fasta = os.path.join(args.directory, args.input_fasta)
    output_dir = args.directory

    # Create nested directories
    try:
        os.makedirs(output_dir)
        print(f"directory '{output_dir}' created successfully.")
    except FileExistsError:
        print(f"directory '{output_dir}' already exist.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{output_dir}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # run pipeline
    seq1_fasta,seq2_fasta = find_and_split_by_pattern(input_fasta, args.pattern, output_dir)

    print(f"Pattern find and split complete. Outputs saved in {output_dir}")


if __name__ == "__main__":
    main()


# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/pattern_find_split.py \
#             -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/aligned_trimmed_sequences.fasta \
#             --pattern [GA].{4}GK[TS] \
#             -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/


