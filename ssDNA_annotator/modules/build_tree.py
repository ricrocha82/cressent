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
        
def sanitize_sequence_names(input_fasta, sanitized_fasta, name_table_file):
    """
    Replace spaces with underscores in sequence names, ensure unique IDs, 
    save the sanitized FASTA file, and create a table with original and sanitized sequence names.
    """
    name_table = []  # To store the old and new names
    sanitized_records = []
    id_map = {}  # To ensure unique IDs

    for record in SeqIO.parse(input_fasta, "fasta"):
        original_id = record.description
        # Replace non-alphanumeric characters with underscores
        sanitized_id = re.sub(r'[^a-zA-Z0-9]', '_', original_id)
        sanitized_id = re.sub('__', '_', sanitized_id)
        
        # Handle duplicate sanitized IDs
        if sanitized_id in id_map:
            id_map[sanitized_id] += 1
            sanitized_id = f"{sanitized_id}_{id_map[sanitized_id]}"
        else:
            id_map[sanitized_id] = 1

        record.id = sanitized_id
        record.description = ""  # Remove description 

        # Add to the name table
        name_table.append({"Original Name": original_id, "Sanitized Name": sanitized_id})
        sanitized_records.append(record)

    # Write the sanitized FASTA file
    SeqIO.write(sanitized_records, sanitized_fasta, "fasta")
    print(f"Sanitized sequence names and saved to {sanitized_fasta}.")

    # Save the name table as a tab-delimited file
    name_table_df = pd.DataFrame(name_table)
    name_table_df.to_csv(name_table_file, sep="\t", index=False)
    print(f"Name table saved to {name_table_file}.")

# def convert_to_phylip(sanitized_fasta, phylip_file):
#     """
#     Convert the sanitized FASTA file to PHYLIP format.
#     """
#     alignment = AlignIO.read(sanitized_fasta, "fasta")
#     AlignIO.write(alignment, phylip_file, "phylip")
#     print(f"Converted {sanitized_fasta} to {phylip_file} in PHYLIP format.")

def generate_phylogenetic_tree(input_fasta):
    """
    Build phylonetic tree using IQTREE
    Step 1: find the appropriate evolutionary model using ModelFinder (Kalyaanamoorthy et al., 2017) > version 1.5.4 and older, -m MFP is the default behavior.
    Step 2: generate tree using bootstrap support generated for 1,000 iterations
    """

    cmd = f"iqtree2 -s {input_fasta} -B 1000 -T AUTO"
    run_command(cmd, "Error building Phylogenetic tree")


def main():
    parser = argparse.ArgumentParser(description="Pipeline for sequence alignment, trimming, motif splitting, and phylogenetic tree generation.")
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing input FASTA file and for saving outputs.")

    args = parser.parse_args()
    prefix = os.path.splitext(os.path.basename(args.input_fasta))[0]
    input_fasta = os.path.join(args.directory, args.input_fasta)
    output_dir = args.directory
    name_table_file = os.path.join(output_dir, f"{prefix}_sanitized_name_table.tsv")
    sanitized_fasta = os.path.join(output_dir, f"{prefix}_sanitized_sequences.fasta")

    # Create nested directories
    try:
        os.makedirs(output_dir)
        print(f"directory '{output_dir}' created successfully.")
    except FileExistsError:
        print(f"directory '{output_dir}' already exist.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{nested_directory}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # run pipeline
    sanitize_sequence_names(input_fasta, sanitized_fasta, name_table_file)
    generate_phylogenetic_tree(sanitized_fasta)

    print(f"Phylogenetic tree analysis complete. Outputs saved in {output_dir}")


if __name__ == "__main__":
    main()


# full Rep sequences
# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
#         -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/aligned_trimmed_sequences.fasta \
#         -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/tree


# Endonuclease domains
# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
#         -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/split_sequences_endonuclease.fasta \
#         -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/tree_endonuc


# Helicase domains
# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
#         -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/split_sequences_helicase.fasta \
#         -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/tree_helic