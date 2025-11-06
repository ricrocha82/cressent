#!/usr/bin/env python3

import argparse
import subprocess
import os
import re
import pandas as pd
from Bio import SeqIO
import logging
import sys

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
        
def sanitize_sequence_names(input_fasta, sanitized_fasta, name_table_file, keep_names=False):
    """
    Process FASTA sequence names:
    - If keep_names=True: Keep only the first word of sequence names
    - If keep_names=False: Replace spaces with underscores
    - In both cases: Ensure unique IDs and create a table with original and sanitized names
    """
    name_table = []  # To store the old and new names
    sanitized_records = []
    id_map = {}  # To ensure unique IDs

    for record in SeqIO.parse(input_fasta, "fasta"):
        original_id = record.description
        
        # Process ID based on parameters
        if keep_names:
            # Keep only the first word before a space
            sanitized_id = original_id.split()[0] if ' ' in original_id else original_id
        else:
            # Replace spaces with underscores
            # sanitized_id = original_id.replace(" ", "_")
            # First replace all non-alphanumeric characters with underscores
            temp_string = re.sub(r'[^a-zA-Z0-9_]', '_', original_id)
            # Then replace multiple consecutive underscores with a single underscore
            sanitized_id = re.sub(r'_+', '_', temp_string)
            # Remove leading/trailing underscores
            sanitized_id = sanitized_id.strip('_')
        
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

    if extra_args:
        # Simply join all arguments with a space
        extra_args_str = " ".join(extra_args)
    else:
        extra_args_str = ""

    
    cmd = f"iqtree -s {input_fasta} -m {model} -B {bootstrap} -T {threads} --prefix {output_prefix} {extra_args_str}"

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
    parser.add_argument("-o", "--output", required=True, default = ".",help="Path to the output directory (Default: working directory)")
    parser.add_argument("-B", "--bootstrap", default=1000, type=int, help="Number of bootstrap iterations (default: 1000)")
    parser.add_argument("-t", "--threads", default="AUTO", help="Number of threads to use (default: AUTO)")
    parser.add_argument("-m", "--model", default="MFP", help="Substitution models (default: MFP - ModelFinder)")
    parser.add_argument("--extra_args", nargs="+", help="Extra arguments to pass directly to IQ-TREE (e.g. --extra_args '-cmax 15')")
    parser.add_argument("--keep_names", action="store_true", help="Keep only first word of sequence IDs")
    

    args = parser.parse_args()
    output_dir = args.output

    prefix = os.path.join(
        output_dir,
        os.path.splitext(
            os.path.basename(args.input_fasta)
            )[0]
    )
    # Determine the input FASTA full path
    def validate_fasta(filename):
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                print("FASTA checked.")
                input_fasta = filename
                return input_fasta
            else:
                sys.exit("Error: Input file is empty or is not in the FASTA format.\n")
                logging.info(f"Using: {input_fasta} is not in the FASTA format")
    # check fasta
    input_fasta = validate_fasta(args.input_fasta)
    log_file = os.path.join(output_dir, "build_tree.log")
    name_table_file = f"{prefix}_sanitized_name_table.tsv"
    sanitized_fasta = f"{prefix}_sanitized_sequences.fasta"
    bootstrap = args.bootstrap
    
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
    # Process sequence names
    sanitize_sequence_names(
        input_fasta, 
        sanitized_fasta, 
        name_table_file,
        keep_names=args.keep_names
    )
    generate_phylogenetic_tree(sanitized_fasta, args.bootstrap, args.threads, prefix, args.extra_args, args.model)
    print(f"Phylogenetic tree analysis complete. Outputs saved in {output_dir}")

    # Log the chosen model
    iqtree_log_file = f"{prefix}.log"
    # iqtree_log_file = output_dir / f"{prefix}.iqtree"
    if args.model == "MFP":
        chosen_model = extract_chosen_model(iqtree_log_file)
    else:
        chosen_model = args.model
    logging.info(f"Chosen evolutionary model for {prefix}: {chosen_model}")


if __name__ == "__main__":
    main()


