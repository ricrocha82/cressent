#!/usr/bin/env python3

import os
import argparse
import subprocess
import logging
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import sys
import re

# Set up logging
def setup_logging(log_file):
    """Configure logging to both file and console."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

# Run shell commands safely
def run_command(command, error_message):
    """Run a shell command and handle errors."""
    logging.info(f"Running command: {command}")
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            text=True, 
            capture_output=True
        )
        if result.stdout:
            logging.info(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
        if e.stderr:
            logging.error(f"Error output: {e.stderr}")
        raise RuntimeError(f"{error_message}: {e}")

def preprocess_fasta(input_fasta, output_fasta, keep_names=False):
    """
    Process FASTA file to ensure unique IDs.
    If keep_names is True, keep only first word of sequence IDs.
    Otherwise, replace spaces with underscores.
    """
    name_table_file = Path(output_fasta).parent / f"{Path(input_fasta).stem}_name_table.tsv"
    logging.info(f"Preprocessing FASTA with keep_names={keep_names}")
    
    name_table = []
    sanitized_records = []
    id_map = {}
    
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        original_id = record.description
        
        # Process ID based on keep_names setting
        if keep_names:
            # Keep only first word before space
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
        
        # Handle duplicate IDs
        if sanitized_id in id_map:
            id_map[sanitized_id] += 1
            sanitized_id = f"{sanitized_id}_{id_map[sanitized_id]}"
        else:
            id_map[sanitized_id] = 1
        
        record.id = sanitized_id
        record.description = ""  # Remove description
        
        name_table.append({"Original_Name": original_id, "Sanitized_Name": sanitized_id})
        sanitized_records.append(record)
    
    # Write sanitized FASTA file
    SeqIO.write(sanitized_records, output_fasta, "fasta")
    logging.info(f"Sanitized {len(sanitized_records)} sequences to {output_fasta}")
    
    # Save name mapping table
    pd.DataFrame(name_table).to_csv(name_table_file, sep="\t", index=False)
    logging.info(f"Name table saved to {name_table_file}")
    
    return output_fasta

def detect_sequence_type(fasta_file):
    """Detect if sequences are nucleotide or protein."""
    nucleotide_chars = {"A", "T", "C", "G", "N"}
    nucleotide_count = 0
    total_chars = 0
    
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            total_chars += len(line)
            nucleotide_count += sum(1 for c in line if c in nucleotide_chars)
    
    if total_chars == 0:
        raise ValueError(f"Empty FASTA file: {fasta_file}")
    
    nucleotide_ratio = nucleotide_count / total_chars
    if nucleotide_ratio > 0.90:
        logging.info("Detected nucleotide sequences (blastn will be used)")
        return "nucleotide"
    else:
        logging.info("Detected protein sequences (blastp will be used)")
        return "protein"

def extract_representative_sequences(input_fasta, rep_sequences, output_fasta):
    """Extract representative sequences to a new FASTA file."""
    logging.info(f"Creating representative sequences FASTA file")
    
    # Count total sequences
    total_count = 0
    kept_count = 0
    
    with open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            total_count += 1
            if record.id in rep_sequences:
                SeqIO.write(record, outfile, "fasta")
                kept_count += 1
    
    logging.info(f"Extracted {kept_count} representative sequences out of {total_count}")
    return output_fasta

def cleanup_temp_dir(temp_dir):
    """Remove temporary files."""
    logging.info("Cleaning up temporary files")
    try:
        for file in Path(temp_dir).glob("*"):
            file.unlink()
        Path(temp_dir).rmdir()
        logging.info("Cleanup completed")
    except Exception as e:
        logging.warning(f"Cleanup error: {e}")

def cluster_sequences(input_fasta, output_dir, threads=32, min_ani=95.0, min_tcov=85.0, 
                     min_qcov=0.0, keep_names=False):
    """Main function to cluster sequences."""
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_file = output_dir / "clustering.log"
    setup_logging(log_file)
    
    # Create temp directory
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    try:
        logging.info("Starting sequence clustering pipeline")
        
        # Define file paths
        if keep_names:
            cleaned_fasta = output_dir / f"{Path(input_fasta).name}"
        else:
            cleaned_fasta = output_dir / f"renamed_{Path(input_fasta).name}"
        db_path = temp_dir / "blast_db"
        blast_output = output_dir / "blast_results.tsv"
        ani_output = output_dir / "ani_results.tsv"
        cluster_output = output_dir / "clusters.tsv"
        cluster_fasta = output_dir / "cluster_sequences.fa"
        
        # Find script paths
        modules_dir = Path(__file__).resolve().parent
        anicalc_script = modules_dir / "anicalc.py"
        aniclust_script = modules_dir / "aniclust.py"
        
        # Check if scripts exist
        if not anicalc_script.exists():
            raise FileNotFoundError(f"anicalc.py not found in {modules_dir}")
        if not aniclust_script.exists():
            raise FileNotFoundError(f"aniclust.py not found in {modules_dir}")
        
        # 1. Preprocess input FASTA
        processed_fasta = preprocess_fasta(input_fasta, cleaned_fasta, keep_names)
        
        # 2. Create BLAST database
        seq_type = detect_sequence_type(processed_fasta)
        db_type = "nucl" if seq_type == "nucleotide" else "prot"
        cmd = f"makeblastdb -in {processed_fasta} -dbtype {db_type} -out {db_path}"
        run_command(cmd, "Failed to create BLAST database")
        
        # 3. Run BLAST search
        if seq_type == "nucleotide":
            blast_cmd = (f"blastn -query {processed_fasta} -db {db_path} -out {blast_output} "
                         f"-outfmt '6 std qlen slen' -max_target_seqs 10000 -num_threads {threads}")
        else:
            blast_cmd = (f"blastp -query {processed_fasta} -db {db_path} -out {blast_output} "
                         f"-outfmt '6 std qlen slen' -evalue 0.00001 -num_threads {threads}")
        run_command(blast_cmd, f"{seq_type} BLAST search failed")
        
        # 4. Calculate ANI
        cmd = f"python {anicalc_script} -i {blast_output} -o {ani_output}"
        run_command(cmd, "ANI calculation failed")
        
        # 5. Cluster sequences
        cmd = (f"python {aniclust_script} --fna {processed_fasta} --ani {ani_output} "
               f"--out {cluster_output} --min_ani {min_ani} --min_tcov {min_tcov} --min_qcov {min_qcov}")
        run_command(cmd, "Sequence clustering failed")
        
        # 6. Process cluster output
        df = pd.read_csv(cluster_output, sep='\t', header=None)
        df.columns = ['Representative_Sequence', 'Sequences']
        rep_sequences = set(df['Representative_Sequence'].tolist())
        df.to_csv(cluster_output, sep='\t', index=False)
        
        logging.info(f"Found {len(rep_sequences)} representative sequences")
        
        # 7. Extract representative sequences
        extract_representative_sequences(processed_fasta, rep_sequences, cluster_fasta)
        
        logging.info(f"Clustering completed successfully! Results saved in {cluster_output}")
        
        return cluster_output, cluster_fasta
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        raise
    finally:
        cleanup_temp_dir(temp_dir)

def main():
    """Parse command line arguments and run sequence clustering."""
    parser = argparse.ArgumentParser(
        description="Sequence clustering using BLAST, anicalc, and aniclust.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input_fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for BLAST (Default = 1)")
    parser.add_argument("--min_ani", type=float, default=95.0, help="Minimum average identity for clustering (Default = 95.0)")
    parser.add_argument("--min_tcov", type=float, default=85.0, help="Minimum target coverage (Default = 85.0)")
    parser.add_argument("--min_qcov", type=float, default=0.0, help="Minimum query coverage (Default = 0.0)")
    parser.add_argument("--keep_names", action="store_true", help="Keep only first word of sequence IDs")
    
    args = parser.parse_args()
    
    # Validate input FASTA
    try:
        with open(args.input_fasta, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if not any(fasta):
                sys.exit("Error: Input file is not in the FASTA format")
            print("FASTA checked.")
    except Exception as e:
        sys.exit(f"Error: Cannot read input file: {e}")
    
    # Run clustering
    try:
        cluster_sequences(
            input_fasta=args.input_fasta,
            output_dir=args.output,
            threads=args.threads,
            min_ani=args.min_ani,
            min_tcov=args.min_tcov,
            min_qcov=args.min_qcov,
            keep_names=args.keep_names
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        exit(1)

if __name__ == "__main__":
    main()