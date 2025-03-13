#!/usr/bin/env python3

import argparse
import subprocess
import os
import re
import pandas as pd
from Bio import SeqIO
import logging

def run_command(command, error_message):
    """
    Helper function to run a shell command and handle errors.
    """
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
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
        name_table.append({"Original_Name": original_id, "Sanitized_Name": sanitized_id})
        sanitized_records.append(record)

    # Write the sanitized FASTA file
    SeqIO.write(sanitized_records, sanitized_fasta, "fasta")
    logging.info(f"Sanitized sequence names and saved to {sanitized_fasta}.")

    # Save the name table as a tab-delimited file
    name_table_df = pd.DataFrame(name_table)
    name_table_df.to_csv(name_table_file, sep="\t", index=False)
    logging.info(f"Name table saved to {name_table_file}.")

# def convert_to_phylip(sanitized_fasta, phylip_file):
#     """
#     Convert the sanitized FASTA file to PHYLIP format.
#     """
#     alignment = AlignIO.read(sanitized_fasta, "fasta")
#     AlignIO.write(alignment, phylip_file, "phylip")
#     print(f"Converted {sanitized_fasta} to {phylip_file} in PHYLIP format.")

def generate_phylogenetic_tree(input_fasta, bootstrap, threads, output_prefix, extra_args, model):
    """
    Build phylonetic tree using IQTREE
    Step 1: find the appropriate evolutionary model using ModelFinder (Kalyaanamoorthy et al., 2017) > version 1.5.4 and older, -m MFP is the default behavior.
    Step 2: generate tree using bootstrap support generated for 1,000 iterations
    """

    cmd = f"iqtree2 -s {input_fasta} -m {model} -B {bootstrap} -T {threads} --prefix {output_prefix} {extra_args}"

    logging.info(f"Running phylogenetic tree generation with command: {cmd}")
    run_command(cmd, "Error building Phylogenetic tree")

# Extract chosen model from IQ-TREE log
def extract_chosen_model(iqtree_log_file):
    with open(iqtree_log_file, 'r') as file:
        for line in file:
            match = re.search(r"Best-fit model according to .+: ([^\s]+)", line)
            if match:
                return match.group(1)
    return "Model not found"


def main():
    parser = argparse.ArgumentParser(description="Pipeline for sequence alignment, trimming, motif splitting, and phylogenetic tree generation.")
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing input FASTA file and for saving outputs.")
    parser.add_argument("-B", "--bootstrap", default=1000, type=int, help="Number of bootstrap iterations (default: 1000)")
    parser.add_argument("-T", "--threads", default="AUTO", help="Number of threads to use (default: AUTO)")
    parser.add_argument("-m", "--model", default="MFP", help="Substitution models (default: MFP - ModelFinder)")
    parser.add_argument("--extra_args", default="", help="Extra arguments to pass directly to IQ-TREE")
    

    args = parser.parse_args()
    prefix = os.path.splitext(os.path.basename(args.input_fasta))[0]
    input_fasta = os.path.join(args.directory, args.input_fasta)
    output_dir = args.directory
    log_file = os.path.join(output_dir, "build_tree.log")
    name_table_file = os.path.join(output_dir, f"{prefix}_sanitized_name_table.tsv")
    sanitized_fasta = os.path.join(output_dir, f"{prefix}_sanitized_sequences.fasta")
    bootstrap = args.bootstrap
    
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
    # Setup logging to save logs in the output directory
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Clear any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    # File handler writes log messages to the log file
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    # Stream handler outputs log messages to the console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logging.info("Starting pipeline")
    sanitize_sequence_names(input_fasta, sanitized_fasta, name_table_file)
    generate_phylogenetic_tree(sanitized_fasta, args.bootstrap, args.threads, prefix, args.extra_args, args.model)
    print(f"Phylogenetic tree analysis complete. Outputs saved in {output_dir}")

    # Log the chosen model
    iqtree_log_file = output_dir / f"{prefix}.iqtree"
    chosen_model = extract_chosen_model(iqtree_log_file)
    logging.info(f"Chosen evolutionary model for {prefix}: {chosen_model}")


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


# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
#         -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/sub_reps_aligned_trimmed_sequences.fasta \
#         -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data

# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_tree.py \
#         -i /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/aligned_db/caps/Microviridae_aligned_cap_trim.fa \
#         -d /fs/scratch/Sullivan_Lab/Ricardo/tree_test