#!/usr/bin/env python

import sys
import os
import csv
import subprocess
import argparse
import logging
import pandas as pd
from Bio import SeqIO, motifs
from Bio.ExPASy import ScanProsite
from Bio import ExPASy
import time


def determine_seq_type(fasta_file):
    """
    Reads the first record in the FASTA file and determines if the sequence is nucleotide or protein.
    """
    record = next(SeqIO.parse(fasta_file, "fasta"))
    nucleotides = set("ATGCUatgcu")
    if set(str(record.seq)).issubset(nucleotides):
        return "dna"
    else:
        return "protein"

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

def run_meme(fasta_file, output_dir, seq_type, meme_args):
    """
    Constructs and executes the MEME command.
    
    Parameters:
    - fasta_file: path to input FASTA file.
    - output_dir: directory where MEME output and additional files will be stored.
    - seq_type: either 'dna' or 'protein'; used to set the proper MEME flag.
    - meme_args: list of additional arguments to pass to MEME.
    """
    # Build the base command
    cmd = ["meme", fasta_file, "-oc", output_dir]
    
    # Add the sequence type flag
    if seq_type == "dna":
        cmd.append("-dna")
    else:
        cmd.append("-protein")
    
    # Add any additional arguments if provided
    if meme_args:
        cmd.extend(meme_args)
    
    # Convert list command into a string for shell execution
    command_str = " ".join(cmd)
    logging.info(f"Running command: {command_str}")
    run_command(command_str, "MEME command failed")

def parse_meme_output(meme_output_dir):
    """
    Parses the MEME output file to extract motif details.
    
    Returns a list of motif dictionaries, e.g.:
     [{'name': 'Motif1', 'consensus': 'ATGC', 'length': 4, 'count': 10, 'instances': [ ... ]}, ...]
    """
    # Parse MEME output using Biopython's Bio.motifs module.
    meme_xml = os.path.join(meme_output_dir, "meme.xml")
    with open(meme_xml) as handle:
        discovered_motifs = motifs.parse(handle, "meme")
    logging.info(f"Parsed {len(list(discovered_motifs))} motifs from {meme_xml}")

    return discovered_motifs

def consensus_motif(motifs, output_dir):
    """
    Extracts a table of consensus sequences from discovered motifs and saves it as a CSV file.
    
    The CSV file includes the following columns:
        - Consensus
        - Length
        - Occurrences
    """
    csv_file = os.path.join(output_dir, "consensus_table.csv")
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        # Write header row
        writer.writerow(["id","consensus", "length", "occurrences"])
        # Write motif data
        for motif in motifs:
            # Assuming each motif object has attributes 'consensus', 'length', and 'num_occurrences'
            writer.writerow([motif.alt_id, motif.consensus, motif.length, motif.num_occurrences])
    logging.info(f"Consensus table saved to: {csv_file}")


def create_motif_table(motifs, output_dir):
    """
    Creates a text file table summarizing motifs for plotting seqlogos.
    The table includes motif name, consensus, length, and count.
    """
    # List to store dictionaries for each instance's attributes
    rows = []

    # Loop over each motif and its instances
    for motif in motifs:
        for instance in motif.instances:
            # Get attribute names: exclude callable attributes and private ones (starting with '_')
            attr_names = [attr for attr in dir(instance) 
                        if not callable(getattr(instance, attr)) and not attr.startswith("_")]
            # Create a dictionary for this instance
            instance_data = {}
            for attr in attr_names:
                try:
                    instance_data[attr] = getattr(instance, attr)
                except Exception as e:
                    instance_data[attr] = None  # In case of any issues accessing an attribute
            rows.append(instance_data)

    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(rows)

    
    # If both 'start' and 'length' columns exist, calculate the 'end' column
    if 'start' in df.columns and 'length' in df.columns:
        df['end'] = df['start'] + df['length']
        # Reorder columns: insert 'end' immediately after 'start'
        cols = list(df.columns)
        if 'end' in cols:
            start_idx = cols.index('start')
            cols.remove('end')
            cols.insert(start_idx + 1, 'end')
            df = df[cols]

    df = df.rename(columns={'motif_name': 'matched'})
    df = df.rename(columns={'sequence_name': 'seqID'})
    motif_table_path = os.path.join(output_dir, "motif_table.csv")
    df.to_csv(motif_table_path, sep='\t', index=False)

    logging.info(f"Motif table saved to: {motif_table_path}")

def create_pwm_matrix(motifs, output_dir):
    """
    For each discovered motif, creates a CSV file containing its PWM matrix.
    The file is named as "<motif.name>_pwm_matrix.csv".
    """
    for motif in motifs:
        # Log the discovered motif details.
        logging.info("Saving PWM matrix")
        
        # Create the file name for this motif's PWM matrix.
        matrix_file = os.path.join(output_dir, f"{motif.alt_id}_pwm_matrix.csv")
        with open(matrix_file, "w", newline="") as f:
            writer = csv.writer(f)
            # Assume pwm is a dictionary: keys are letters and values are lists of numbers.
            # Determine the number of positions from the length of the first entry.
            positions = list(range(1, len(next(iter(motif.pwm.values()))) + 1))
            # Write header row: an empty cell then columns for each position.
            header = [""] + [f"Pos {i}" for i in positions]
            writer.writerow(header)
            # Write each letter's PWM row.
            for letter in sorted(motif.pwm.keys()):
                row = [letter] + list(motif.pwm[letter])
                writer.writerow(row)
        logging.info(f"PWM matrix for motif {motif.name} saved to: {matrix_file}")

def run_scanprosite(fasta_file, seq_type):
    """
    Scans each protein sequence in the input FASTA file using Biopython's ScanProsite.
    For each record, a request is sent to the ExPASy server and the scan results are parsed.
    A DataFrame is built for each sequence and then combined into a single DataFrame,
    which is saved as a tab-separated CSV file in the output directory.
    
    Parameters:
      fasta_file: Path to the input FASTA file containing protein sequences.
      output_dir: Directory where the scan results CSV file will be saved.
    
    Returns:
      combined_df: The combined pandas DataFrame with scan results for all sequences.
    """
    # Add the sequence type flag.
    if seq_type == "dna":
        logging.info("ScanProsite requires protein sequences. Skipping ScanProsite analysis.")
        return pd.DataFrame()  # Return an empty DataFrame for DNA sequences.
    else:
        all_dfs = []
        logging.info("Scanning protein sequences online against the Prosite database")

        for record in SeqIO.parse(fasta_file, "fasta"):
            # Prepare the FASTA string for the protein sequence.
            fasta_str = record.format("fasta")
            try:
                # Call ScanProsite via Biopython; this sends a request to the ExPASy server.
                handle = ScanProsite.scan(seq=fasta_str, output = 'xml', skip=0)
                scan_results = ScanProsite.read(handle)
                # Build a DataFrame from the scan results.
                df = pd.DataFrame(scan_results)
                # Add columns to track which sequence produced the results.
                df["record_id"] = record.id
                df["sequence_name"] = record.description
                all_dfs.append(df)
            except Exception as e:
                logging.error(f"Error scanning sequence {record.id}: {e}")
            # Pause briefly to avoid overwhelming the ExPASy server.
            time.sleep(1)

        # Combine all individual DataFrames into one.
        if all_dfs:
            combined_df = pd.concat(all_dfs, ignore_index=True)
        else:
            combined_df = pd.DataFrame()

        return combined_df

def add_annotation_to_scanprosite_results(df, output_dir):
    """
    Adds an 'annotation' column to the ScanProsite results DataFrame.
    
    For each unique signature_ac in the DataFrame, this function:
      1. Uses Bio.ExPASy's get_prosite_raw() to fetch the raw Prosite entry.
      2. Parses the raw text to extract the description line (the line starting with "DE").
      3. Stores the description as the annotation for that signature.
      
    Parameters:
      df (pandas.DataFrame): DataFrame containing the ScanProsite results. 
                             It must include a 'signature_ac' column.
                             
    Returns:
      pandas.DataFrame: The input DataFrame with an added 'annotation' column.
    """
    annotations = {}

    logging.info("Adding Prosite annotation")

    if df.empty:
        logging.info("ScanProsite results DataFrame is empty. Skipping annotation step.")
        return df
    
    # Iterate over each unique signature_ac in the DataFrame.
    for signature in df["signature_ac"].unique():
        try:
            # Fetch the raw Prosite entry for this signature.
            handle = ExPASy.get_prosite_raw(signature)
            text = handle.read()
            # Initialize annotation as an empty string.
            annotation = ""
            # Parse the raw text: the description is on the line starting with "DE".
            for line in text.splitlines():
                if line.startswith("DE"):
                    # Extract the description (remove the "DE" prefix and any trailing semicolon).
                    annotation = line[2:].strip().rstrip(";")
                    break
            annotations[signature] = annotation
        except Exception as e:
            annotations[signature] = ""
            logging.error(f"Error retrieving annotation for {signature}: {e}")
    
    # Map the annotations to a new 'annotation' column.
    df["prosite_ann"] = df["signature_ac"].map(annotations)

    try:
        # Save the updated DataFrame with annotations.
        output_csv = os.path.join(output_dir, "scanprosite_results.csv")
        df.to_csv(output_csv, sep='\t', index=False)
        logging.info(f"ScanProsite results saved to: {output_csv}")
    except Exception as e:
        logging.error(f"Error saving annotated ScanProsite results: {e}")


def main():
    parser = argparse.ArgumentParser(description="motif_discovery: Discover de novo motifs using MEME")
    parser.add_argument("-i","--input_fasta", help="Input FASTA file", required=True)
    parser.add_argument("-o", "--output", default=".", help="Path to the output directory (Default: working directory) (used for MEME and generated files)")
    parser.add_argument("-nmotifs", type=int, default=3, help="Number of motifs to find (Default = 3)")
    parser.add_argument("-minw", type=int, default=5, help="Minimum motif width (Default = 5)")
    parser.add_argument("-maxw", type=int, default=10, help="Maximum motif width (Default = 10)")
    parser.add_argument("--meme_extra", nargs="+", help="Additional MEME arguments (list format)")
    parser.add_argument("--scanprosite", action='store_true', help='Run ScanProsite')
    
    args = parser.parse_args()
    # Check if the output directory exists.
    if os.path.exists(args.output):
        print(f"Warning: Output directory '{args.output}' already exists. Its contents will be overwritten.")
    else:
        os.makedirs(args.output, exist_ok=True)
    
    
    # Configure logging to save logs in the output directory.
    log_file = os.path.join(args.output, "motif_discovery.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    logging.info("Starting motif_discovery module.")
    
    seq_type = determine_seq_type(args.input_fasta)
    logging.info(f"Detected sequence type: {seq_type}")
    
    # Set default MEME arguments
    default_args = ["-nmotifs", str(args.nmotifs), "-minw", str(args.minw), "-maxw", str(args.maxw)]
    
    # Use extra arguments if provided, otherwise use defaults
    if args.meme_extra:
        meme_args = default_args + args.meme_extra
        logging.info(f"Using provided MEME arguments: {meme_args}")
    else:
        meme_args = default_args
        logging.info(f"Using default MEME arguments: {meme_args}")
    
    run_meme(args.input_fasta, args.output, seq_type, meme_args)
    
    # Parse MEME output (assumes MEME XML output is stored in <output>/meme.xml).
    motifs = parse_meme_output(args.output)
    
    create_motif_table(motifs, args.output)
    create_pwm_matrix(motifs, args.output)
    consensus_motif(motifs, args.output)


    # Determine if ScanProsite will be used
    scanprosite = args.scanprosite
    if scanprosite:
        logging.info("ScanProsite selected.")
        # Run ScanProsite and get the combined DataFrame.
        scan_df = run_scanprosite(args.input_fasta, seq_type)
        # Add the annotation column.
        add_annotation_to_scanprosite_results(scan_df, args.output)


if __name__ == "__main__":
    main()


# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/motif_disc.py \
#                     -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/seq_meme.fa \
#                     -o /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc
