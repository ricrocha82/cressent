#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
import logging
import pandas as pd
import re
from Bio import SeqIO
from typing import List, Optional


def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "alignment.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def run_command(command: str, error_message: str) -> None:
    """Helper function to run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
        exit(1)

def get_database_files(db_path: str, families: List[str], protein_type: Optional[str] = None) -> List[str]:
    """Get relevant database files based on family names and protein type."""
    db_files = []
    
    if "all" in families:  
        families = [f.split(".")[0] for f in os.listdir(db_path) if f.endswith(".fa") and f != "all.fa"]

    for family in families:
        if protein_type:
            filename = f"{protein_type}/{family}.fa"
        else:
            filename = f"{family}.fa"
        
        file_path = os.path.join(db_path, filename)
        if os.path.exists(file_path):
            db_files.append(file_path)
        else:
            logging.warning(f"Database file not found: {file_path}")

    return db_files

def merge_fasta_files(input_file: str, db_files: List[str], output_file: str) -> None:
    """Merge input FASTA file with database files."""
    with open(output_file, 'w') as out_f:
        # First, copy input file
        if os.path.exists(input_file):
            with open(input_file, 'r') as in_f:
                shutil.copyfileobj(in_f, out_f)
            logging.info(f"Added sequences from {input_file}")
        
        # Then add database sequences
        for db_file in db_files:
            if os.path.exists(db_file):
                with open(db_file, 'r') as in_f:
                    shutil.copyfileobj(in_f, out_f)
                logging.info(f"Added sequences from {db_file}")

def create_metadata_dataframe(fasta_file: str, family: str = "input", source: str = "input") -> pd.DataFrame:
    """Create a metadata dataframe from a FASTA file with debug logging."""
    metadata = []
    logging.info(f"Starting to process metadata for file: {fasta_file}")
    
    try:
        # Check if file exists
        if not os.path.exists(fasta_file):
            logging.error(f"FASTA file not found: {fasta_file}")
            return pd.DataFrame()
        
        # Try to parse sequences
        sequence_count = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_count += 1
            protein_id = record.id
            protein_description = record.description
            
            logging.debug(f"Processing sequence: {protein_id}")
            
            # Extract scientific name
            scientific_name_match = re.search(r'\[(.*?)\]', protein_description)
            scientific_name = scientific_name_match.group(1) if scientific_name_match else "Unknown"
            
            # Extract protein name
            protein_name = re.sub(r' \[.*\]', '', protein_description).split(' ', 1)[-1]
            
            metadata.append({
                "protein_id": protein_id,
                "protein_description": protein_description,
                "family": family,
                "scientific_name": scientific_name,
                "protein_name": protein_name,
                "source": source
            })
            
            logging.debug(f"Added metadata for: {protein_id}")
        
        logging.info(f"Processed {sequence_count} sequences from {fasta_file}")
        
        if not metadata:
            logging.warning(f"No sequences were processed from {fasta_file}")
            return pd.DataFrame()
            
        df = pd.DataFrame(metadata)
        logging.info(f"Created DataFrame with {len(df)} rows")
        return df
        
    except Exception as e:
        logging.error(f"Error processing {fasta_file}: {str(e)}")
        return pd.DataFrame()

def save_metadata(input_fasta: str, db_files: List[str], output_dir: str) -> None:
    """Generate and save metadata with enhanced error checking and logging."""
    logging.info("Starting metadata creation process")
    metadata_dfs = []
    
    # Process input sequences
    logging.info(f"Processing input FASTA file: {input_fasta}")
    input_df = create_metadata_dataframe(input_fasta, family="my_samples", source="input")
    if not input_df.empty:
        metadata_dfs.append(input_df)
        logging.info(f"Successfully created metadata for input file with {len(input_df)} entries")
    else:
        logging.warning("No metadata created for input file")
    
    # Process database sequences if any
    if db_files:
        for db_file in db_files:
            logging.info(f"Processing database file: {db_file}")
            family_name = os.path.splitext(os.path.basename(db_file))[0]
            db_df = create_metadata_dataframe(db_file, family=family_name, source="database")
            if not db_df.empty:
                metadata_dfs.append(db_df)
                logging.info(f"Successfully created metadata for {db_file} with {len(db_df)} entries")
    
    # Save combined metadata
    if metadata_dfs:
        try:
            combined_df = pd.concat(metadata_dfs, ignore_index=True)
            metadata_file = os.path.join(output_dir, "metadata.csv")
            combined_df.to_csv(metadata_file, index=False)
            logging.info(f"Saved metadata to {metadata_file} with {len(combined_df)} total entries")
        except Exception as e:
            logging.error(f"Error saving metadata: {str(e)}")
    else:
        logging.error("No metadata was created - no valid data found in any input files")

def align_sequences(input_fasta: str, output_dir: str, prefix: str, mafft_ep: float = 0.123, threads: int = 1) -> str:
    """Align sequences using MAFFT and save output in the output directory."""
    logging.info("Starting alignment using MAFFT")
    aligned_file = os.path.join(output_dir, f"{prefix}_aligned_sequences.fasta")
    cmd = f"mafft --localpair --maxiterate 1000 --leavegappyregion --thread {threads} --ep {mafft_ep} {input_fasta} > {aligned_file}"
    run_command(cmd, "Error running MAFFT")
    logging.info(f"Alignment saved to {aligned_file}")
    return aligned_file

def trim_alignment(aligned_fasta: str, output_dir: str, prefix: str, gap_threshold: float = 0.2) -> str:
    """Trim alignment using TrimAl."""
    logging.info("Starting trimming using TrimAl")
    trimmed_file = os.path.join(output_dir, f"{prefix}_aligned_trimmed_sequences.fasta")
    cmd = f"trimal -keepheader -in {aligned_fasta} -out {trimmed_file} -gt {gap_threshold}"
    run_command(cmd, "Error running TrimAl")
    logging.info(f"Trimmed alignment saved to {trimmed_file}")
    return trimmed_file

def check_dependency(command: str) -> bool:
    """Check if a command-line tool is available."""
    try:
        subprocess.run([command, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        logging.info(f"{command} is installed and available.")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error(f"Error: {command} is not installed or not in PATH.")
        return False

def main():
    parser = argparse.ArgumentParser(description="Pipeline for sequence alignment and trimming.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("-i", "--input_fasta", required=True, help="Input FASTA file with sequences")
    parser.add_argument("-d", "--directory", default=".", help="Directory for saving outputs")
    parser.add_argument("--mafft_ep", type=float, default=0.123, help="Alignment length for MAFFT (default: 0.123)")
    parser.add_argument("--gap_threshold", type=float, default=0.2, help="Gap threshold for TrimAl (default: 0.2)")
    parser.add_argument("--db_family", nargs='+', help="List of family names or 'all' to use multiple database sequences")
    parser.add_argument("--db_path", default="./db", help="Path to the database FASTA files")
    parser.add_argument("--protein_type", choices=['reps', 'caps', 'orf'], help="Specify protein type (Rep or Cap) for database files")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.directory, exist_ok=True)
    setup_logging(args.directory)
    logging.info(f"Output directory set to {args.directory}")

    # validate dependencies
    if not check_dependency("mafft") or not check_dependency("trimal"):
        logging.error("Required dependencies (MAFFT, TrimAl) are missing. Exiting.")
        exit(1)

    
    # Set up input and output paths
    prefix = os.path.splitext(os.path.basename(args.input_fasta))[0]
    input_fasta = args.input_fasta

    # Generate metadata only for the input file first
    save_metadata(args.input_fasta, [], args.directory)  # No database sequences yet
    
     # Handle database integration if requested
    db_files = []
    # Handle database integration if requested
    if args.db_family:
        db_files = get_database_files(args.db_path, args.db_family, args.protein_type)
        
        # Generate metadata for database sequences separately
        save_metadata(args.input_fasta, db_files, args.directory)
        
        # Merge input and database sequences
        merged_fasta = os.path.join(args.directory, f"{prefix}_merged.fasta")
        merge_fasta_files(args.input_fasta, db_files, merged_fasta)
        input_fasta = merged_fasta  # Use merged file for alignment
    else:
        input_fasta = args.input_fasta  # No database, just use input

    # Run alignment and trimming after save the metada
    aligned_fasta = align_sequences(input_fasta, args.directory, prefix, 
                                        mafft_ep=args.mafft_ep, threads=args.threads)
    trim_alignment(aligned_fasta, args.directory, prefix, 
                    gap_threshold=args.gap_threshold)

if __name__ == "__main__":
    main()


#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/

#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_test/test_data/Circoviridae_caps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_3/

#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_test/test_data/Circoviridae_reps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_3/

#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_test/test_data/Circoviridae_reps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_3/


#/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/concat/reps


# with DB
# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py  --threads 24 \
#                     --input_fasta /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/sub_reps.fa \
#                     --db_family Circoviridae \
#                     --db_path /fs/project/PAS1117/ricardo/ssDNA_tool/DB \
#                      --protein_type reps \
#                     -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/align_test

# only input sequences
# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py  --threads 24 \
#                       --input_fasta /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/sub_reps.fa \
#                   -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_4/