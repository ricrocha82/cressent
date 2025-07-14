import os
import re
import argparse
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def setup_logging(output_dir: str) -> None:
    """Setup logging configuration with log file in output directory"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    log_file = os.path.join(output_dir, 'adjust_seq.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logging.info(f"Log file created at: {log_file}")

def detect_sequence_type(sequence: str) -> str:
    """Detect if a sequence is nucleotide or protein based on its composition."""
    sequence_upper = sequence.upper()
    
    # Check for nucleotide sequence (A, T, G, C, N and possibly U for RNA)
    nucleotide_chars = set('ATGCNU')
    protein_chars = set('ACDEFGHIKLMNPQRSTVWY')
    
    # Count valid nucleotide and protein characters
    nucleotide_count = sum(1 for char in sequence_upper if char in nucleotide_chars)
    protein_count = sum(1 for char in sequence_upper if char in protein_chars)
    
    # Calculate percentages
    total_length = len(sequence_upper)
    nucleotide_percentage = nucleotide_count / total_length
    protein_percentage = protein_count / total_length
    
    # Decision logic
    if nucleotide_percentage >= 0.95:  # 95% or more nucleotide characters
        return "nucleotide"
    elif protein_percentage >= 0.90:  # 90% or more protein characters
        return "protein"
    else:
        return "unknown"

def adjust_sequence_start(sequence: str, motif: str, seq_type: str) -> str:
    """Rotate the sequence so it starts with the given motif."""
    
    # Use regex search to find the motif
    match = re.search(motif, sequence)
    if match:
        motif_index = match.start()
        logging.info(f"Motif '{motif}' found at position {motif_index} in {seq_type} sequence. Adjusting sequence start position.")
        return sequence[motif_index:] + sequence[:motif_index]
    else:
        logging.warning(f"Motif '{motif}' not found in {seq_type} sequence: {sequence[:30]}...")
        return None  # Return None if motif not found

def process_fasta(input_fasta: str, output_dir: str, motif: str) -> str:
    """Adjust sequences in a FASTA file to start with the specified motif."""
    
    # Generate output file path
    prefix = os.path.splitext(os.path.basename(input_fasta))[0]
    output_fasta = os.path.join(output_dir, f"{prefix}_motif_adj.fa")
    
    logging.info(f"Processing FASTA file: {input_fasta}")
    logging.info(f"Output file: {output_fasta}")
    logging.info(f"Motif pattern: {motif}")
    
    sequences_processed = 0
    sequences_adjusted = 0
    sequences_skipped = 0
    
    nucleotide_count = 0
    protein_count = 0
    unknown_count = 0
    
    with open(output_fasta, "w") as out_fh:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_str = str(record.seq)
            sequences_processed += 1
            
            # Detect sequence type
            seq_type = detect_sequence_type(seq_str)
            
            if seq_type == "nucleotide":
                nucleotide_count += 1
            elif seq_type == "protein":
                protein_count += 1
            else:
                unknown_count += 1
            
            logging.info(f"Processing sequence {record.id}: {seq_type} sequence (length: {len(seq_str)})")
            
            # Adjust sequence start
            adjusted_seq = adjust_sequence_start(seq_str, motif, seq_type)
            
            if adjusted_seq is not None:
                record.seq = Seq(adjusted_seq)
                SeqIO.write(record, out_fh, "fasta")
                sequences_adjusted += 1
            else:
                sequences_skipped += 1
                logging.warning(f"Skipping sequence {record.id}: motif not found")
    
    # Print summary
    logging.info("Processing Summary:")
    logging.info(f"Total sequences processed: {sequences_processed}")
    logging.info(f"Sequences adjusted: {sequences_adjusted}")
    logging.info(f"Sequences skipped: {sequences_skipped}")
    logging.info(f"Nucleotide sequences: {nucleotide_count}")
    logging.info(f"Protein sequences: {protein_count}")
    logging.info(f"Unknown sequences: {unknown_count}")
    
    return output_fasta

def validate_fasta(filename: str) -> str:
    """Validate that the input file is in FASTA format."""
    logging.info(f"Validating FASTA file: {filename}")
    
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                logging.info("FASTA file validated successfully.")
                return filename
            else:
                logging.error("Input file is not in valid FASTA format.")
                sys.exit("Error: Input file is not in valid FASTA format.\n")
    except FileNotFoundError:
        logging.error(f"Input file not found: {filename}")
        sys.exit(f"Error: Input file not found: {filename}\n")
    except Exception as e:
        logging.error(f"Error validating FASTA file: {str(e)}")
        sys.exit(f"Error validating FASTA file: {str(e)}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Adjust sequences in a FASTA file to start with a specified motif. "
                   "Supports both nucleotide and protein sequences."
    )
    parser.add_argument(
        "-i", "--input_fasta", 
        required=True,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-o", "--output", 
        default=".", 
        help="Path to the output directory. (Default: current working directory)"
    )
    parser.add_argument(
        "-m", "--motif", 
        default="TAGTATTAC", 
        help="Motif/regex pattern to adjust sequences to start with (default: TAGTATTAC)."
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup logging
    setup_logging(output_dir)
    
    logging.info(f"Input FASTA: {args.input_fasta}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Motif pattern: {args.motif}")
    
    try:
        # Validate input FASTA file
        input_fasta = validate_fasta(args.input_fasta)
        
        # Process the FASTA file
        output_file = process_fasta(input_fasta, output_dir, args.motif)
        
        logging.info(f"Processing completed successfully!")
        logging.info(f"Output file created: {output_file}")
        
    except Exception as e:
        logging.error(f"An error occurred during processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
