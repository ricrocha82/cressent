#!/usr/bin/env python3
"""
CRESS DNA/Protein Virus Decontamination Pipeline

This script screens sequence data against a database of known viral contaminants
and filters out potential contamination. It can handle both nucleotide and protein sequences.

Usage:
    python decontaminate.py --input_fasta sequences.fasta --db viralContaminants.fasta --output-dir results

Dependencies:
    - Biopython
    - pandas
    - BLAST+ command line tools
"""

import argparse
import os
import sys
import logging
import subprocess
import tempfile
from typing import Dict, List, Tuple, Set, Optional
import pandas as pd
from Bio import SeqIO

# Create a global logger object
logger = logging.getLogger('viral_decontamination')

def setup_logging(log_file: str) -> logging.Logger:
    """
    Configure logging to both file and console.
    
    Args:
        log_file: Path to the log file
        
    Returns:
        Configured logger object
    """
    global logger
    
    # Reset logger if it already has handlers
    if logger.handlers:
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
    
    logger.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

def check_dependencies() -> bool:
    """Check if required external dependencies are installed."""
    try:
        # Check for blastn
        subprocess.run(["blastn", "-version"], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE, 
                      check=True)
        
        # Check for blastp
        subprocess.run(["blastp", "-version"], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE, 
                      check=True)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        logger.error("BLAST+ tools not found. Please install BLAST+ and ensure it's in your PATH.")
        return False

def validate_fasta(filename):
    """Validate that the input file is in FASTA format."""
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                logger.info(f"FASTA file validated: {filename}")
                return filename
            else:
                logger.error(f"Error: Input file is not in the FASTA format: {filename}")
                raise ValueError("Input file is not in the FASTA format")
    except Exception as e:
        logger.error(f"Error validating FASTA file: {e}")
        raise

def detect_sequence_type(fasta_file: str) -> str:
    """
    Detect whether the input FASTA file contains nucleotide or protein sequences.
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        'nucl' or 'prot' based on detection
    """
    # Read first few sequences to determine type
    nuc_chars = set('ATGCNatgcn')
    prot_chars = set('ACDEFGHIKLMNPQRSTVWYBXZacdefghiklmnpqrstvwybxz') - nuc_chars # Amino acids not in nucleotide alphabet
    
    try:
        sequences = []
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append(str(record.seq))
                if len(sequences) >= 5:  # Check up to 5 sequences
                    break
        
        # If we couldn't get any sequences, default to nucleotide
        if not sequences:
            logger.warning("Couldn't detect sequence type. Defaulting to nucleotide.")
            return 'nucl'
        
        # Check each sequence for protein-specific amino acids
        is_protein = False
        for seq in sequences:
            seq_chars = set(seq.upper())
            if seq_chars.intersection(prot_chars):
                is_protein = True
                break
        
        if is_protein:
            logger.info("Detected protein sequences")
            return 'prot'
        else:
            logger.info("Detected nucleotide sequences")
            return 'nucl'
            
    except Exception as e:
        logger.error(f"Error detecting sequence type: {e}")
        raise RuntimeError(f"Failed to detect sequence type: {e}")

def run_blast(query_file: str, db_file: str, output_file: str, 
              seq_type: str = 'nucl',
              evalue: float = 1e-10, 
              identity: float = 90.0,
              threads: int = 1) -> bool:
    """
    Run BLASTN or BLASTP to identify potential contaminants.
    
    Args:
        query_file: Path to the input FASTA file
        db_file: Path to the contaminant database FASTA file
        output_file: Path to save BLAST results
        seq_type: Type of sequences ('nucl' or 'prot')
        evalue: E-value threshold (default: 1e-10)
        identity: Percent identity threshold (default: 90.0)
        threads: Number of CPU threads to use
        
    Returns:
        True if BLAST ran successfully, False otherwise
    """
    blast_program = "blastn" if seq_type == 'nucl' else "blastp"
    logger.info(f"Running {blast_program} against contaminant database...")
    
    # Create a temporary BLAST database if needed
    db_path = os.path.splitext(db_file)[0]
    if not os.path.exists(f"{db_path}.{('nhr' if seq_type == 'nucl' else 'phr')}"):
        logger.info(f"Creating BLAST database from {db_file}")
        cmd_makedb = [
            "makeblastdb", 
            "-in", db_file,
            "-dbtype", seq_type,
            "-out", db_path
        ]
        try:
            subprocess.run(cmd_makedb, check=True, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE)
        except subprocess.SubprocessError as e:
            logger.error(f"Failed to create BLAST database: {e}")
            return False
    
    # Run BLAST
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
    cmd_blast = [
        blast_program,
        "-query", query_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", outfmt,
        "-evalue", str(evalue),
        "-num_threads", str(threads)
    ]
    
    # Add percent identity only for blastn (blastp uses different parameters)
    if seq_type == 'nucl':
        cmd_blast.extend(["-perc_identity", str(identity)])
    
    try:
        subprocess.run(cmd_blast, check=True, 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE)
        logger.info(f"BLAST search completed, results in {output_file}")
        return True
    except subprocess.SubprocessError as e:
        logger.error(f"BLAST search failed: {e}")
        return False

def parse_blast_results(blast_output: str, 
                        min_identity: float = 90.0,
                        min_coverage: float = 50.0) -> Set[str]:
    """
    Parse BLAST results to identify contaminant sequences.
    
    Args:
        blast_output: Path to the BLAST output file
        min_identity: Minimum percent identity to consider a match (default: 90.0)
        min_coverage: Minimum query coverage to consider a match (default: 50.0)
        
    Returns:
        Set of sequence IDs identified as contaminants
    """
    # Check if the output file exists and has content
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        logger.info("No BLAST hits found against contaminant database.")
        return set()
    
    # Define column names for BLAST output format 6 with added qcovs
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"
    ]
    
    # Parse BLAST output
    try:
        df = pd.read_csv(blast_output, sep='\t', header=None, names=columns)
        
        # Filter by identity and coverage thresholds
        filtered_df = df[(df['pident'] >= min_identity) & (df['qcovs'] >= min_coverage)]
        
        # Get unique sequence IDs that match contaminants
        contaminant_ids = set(filtered_df['qseqid'].unique())
        
        logger.info(f"Identified {len(contaminant_ids)} sequences as potential contaminants")
        return contaminant_ids
        
    except Exception as e:
        logger.error(f"Error parsing BLAST output: {e}")
        return set()

def filter_sequences(input_file: str, 
                     output_file: str, 
                     contaminant_ids: Set[str],
                     stats_file: Optional[str] = None) -> Tuple[int, int]:
    """
    Filter out contaminant sequences from the input file.
    
    Args:
        input_file: Path to the input FASTA file
        output_file: Path to save filtered sequences
        contaminant_ids: Set of sequence IDs to remove
        stats_file: Optional path to save statistics (default: None)
        
    Returns:
        Tuple of (total_sequences, filtered_sequences)
    """
    logger.info(f"Filtering contaminant sequences...")
    
    # Read and filter sequences
    total_count = 0
    filtered_count = 0
    filtered_sequences = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        total_count += 1
        if record.id not in contaminant_ids:
            filtered_sequences.append(record)
        else:
            filtered_count += 1
    
    # Write filtered sequences
    SeqIO.write(filtered_sequences, output_file, "fasta")
    
    logger.info(f"Removed {filtered_count} contaminant sequences out of {total_count}")
    
    # Write statistics if requested
    if stats_file:
        with open(stats_file, 'w') as f:
            f.write(f"Total sequences: {total_count}\n")
            f.write(f"Identified contaminants: {filtered_count}\n")
            f.write(f"Clean sequences: {total_count - filtered_count}\n")
            f.write(f"Contamination rate: {filtered_count/total_count*100:.2f}%\n")
            
            if filtered_count > 0:
                f.write("\nContaminant sequence IDs:\n")
                for seq_id in sorted(contaminant_ids):
                    f.write(f"{seq_id}\n")
    
    return total_count, filtered_count

def decontaminate(input_file: str, 
                  db_file: str, 
                  output_file: str,
                  seq_type: str = None,
                  evalue: float = 1e-10,
                  identity: float = 90.0, 
                  coverage: float = 50.0,
                  threads: int = 1,
                  keep_temp: bool = False,
                  blast_output: str = None,
                  stats_file: str = None) -> None:
    """
    Main decontamination function.
    
    Args:
        input_file: Path to the input FASTA file
        db_file: Path to the contaminant database
        output_file: Path to save filtered sequences
        seq_type: Type of sequences ('nucl' or 'prot', auto-detect if None)
        evalue: BLAST E-value threshold (default: 1e-10)
        identity: Percent identity threshold (default: 90.0)
        coverage: Query coverage threshold (default: 50.0)
        threads: Number of CPU threads to use (default: 1)
        keep_temp: Whether to keep temporary files (default: False)
        blast_output: Path to save BLAST results (default: None)
        stats_file: Path to save statistics (default: None)
    """
    logger.info(f"Starting decontamination of {input_file}")
    
    # Auto-detect sequence type if not specified
    if seq_type is None:
        seq_type = detect_sequence_type(input_file)
    
    # Run BLAST to identify contaminants
    if not run_blast(input_file, db_file, blast_output, seq_type, evalue, identity, threads):
        logger.error("BLAST search failed. Exiting.")
        sys.exit(1)
    
    # Parse BLAST results
    contaminant_ids = parse_blast_results(blast_output, identity, coverage)
    
    # Filter sequences
    total, filtered = filter_sequences(input_file, output_file, contaminant_ids, stats_file)
    
    logger.info(f"Decontamination complete: {filtered}/{total} sequences identified as contaminants")
    
    # Remove BLAST results if not keeping them
    if not keep_temp and blast_output != stats_file and os.path.exists(blast_output):
        try:
            os.remove(blast_output)
            logger.info(f"Removed temporary BLAST results file")
        except Exception as e:
            logger.warning(f"Failed to remove temporary file {blast_output}: {e}")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Filter viral contaminants from sequence data")
    
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file")
    parser.add_argument("--db", required=True, help="Contaminant database FASTA file")
    parser.add_argument("-o","--output", required=True, help="Path to the output directory (Default: working directory)")
    parser.add_argument("--output-name", default="clean_sequences", 
                        help="Base name for output files (default: clean_sequences)")
    parser.add_argument("--seq-type", choices=['nucl', 'prot'], 
                        help="Sequence type (auto-detect if not specified)")
    parser.add_argument("--evalue", type=float, default=1e-10, 
                        help="BLAST E-value threshold (default: 1e-10)")
    parser.add_argument("--identity", type=float, default=90.0, 
                        help="Minimum percent identity to consider a match (default: 90.0)")
    parser.add_argument("--coverage", type=float, default=50.0, 
                        help="Minimum query coverage to consider a match (default: 50.0)")
    parser.add_argument("-t","--threads", type=int, default=1, 
                        help="Number of CPU threads for BLAST (default: 1)")
    parser.add_argument("--keep-temp", action="store_true", 
                        help="Keep temporary BLAST output files")
    
    return parser.parse_args()

def main():
    """Main function to run the contamination detection pipeline (for CLI import)."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # Determine the input FASTA full path
    input_fasta = validate_fasta(args.input_fasta)
    
    # Set up file paths
    output_fasta = os.path.join(output_dir, f"{args.output_name}")
    stats_file = os.path.join(output_dir, f"{args.output_name}_stats.txt")
    blast_output_file = os.path.join(output_dir, f"{args.output_name}_blast.tsv")
    log_file = os.path.join(output_dir, f"{args.output_name}_decontamination.log")
    
    # Set up logging - make logger global
    global logger
    logger = setup_logging(log_file)
    logger.info(f"Starting viral decontamination process")
    logger.info(f"Input file: {input_fasta}")
    logger.info(f"Contaminant database: {args.db}")
    logger.info(f"Output directory: {output_dir}")
    
    # Check dependencies
    if not check_dependencies():
        logger.error("Required dependencies not found. Exiting.")
        return 1
    
    # Create temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define temporary BLAST output file if not keeping results
        if args.keep_temp:
            blast_result_file = blast_output_file
        else:
            blast_result_file = os.path.join(temp_dir, "contaminant_hits.tsv")
        
        # Run decontamination
        decontaminate(
            input_fasta, 
            args.db, 
            output_fasta,
            args.seq_type,
            args.evalue,
            args.identity,
            args.coverage,
            args.threads,
            args.keep_temp,
            blast_result_file,
            stats_file
        )
    
    logger.info(f"Decontamination process completed")
    logger.info(f"Results saved to {output_dir}")
    logger.info(f"- Filtered sequences: {output_fasta}")
    logger.info(f"- Statistics file: {stats_file}")
    if args.keep_temp:
        logger.info(f"- BLAST results: {blast_output_file}")
    logger.info(f"- Log file: {log_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

