#!/usr/bin/env python3
"""
CRESS DNA Virus Contaminant Database Builder

This script builds a comprehensive database of viral contaminants from a CSV file
of accession numbers for use in decontamination pipelines.

Usage:
    python build_contaminant_db.py --accession-csv decont_accesion_list.csv --output viralContaminants.fasta

Dependencies:
    - Biopython
    - pandas
"""

import argparse
import os
import sys
import logging
import time
from typing import Dict, List, Optional
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Function to set up logging to file and console
def setup_logging(log_file):
    """Configure logging to both file and console."""
    # Create logger
    logger = logging.getLogger('contaminant_db_builder')
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

def set_entrez_email(email: str) -> None:
    """Set email for Entrez to track usage."""
    Entrez.email = email

def download_sequence_by_accession(accession: str, retries: int = 3) -> Optional[SeqRecord]:
    """
    Download a sequence from NCBI using its accession number.
    Includes retry logic for reliability.
    """
    for attempt in range(retries):
        try:
            logger.info(f"Downloading sequence: {accession}")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Add the accession to the description for better tracking
            if not record.description:
                record.description = f"Viral contaminant {accession}"
            else:
                record.description = f"{record.description} [Contaminant]"
                
            return record
        except Exception as e:
            logger.warning(f"Attempt {attempt+1}/{retries} failed for {accession}: {str(e)}")
            if attempt < retries - 1:
                time.sleep(2)  # Wait 2 seconds before retrying
            else:
                logger.error(f"Failed to download {accession} after {retries} attempts")
                return None

def load_accessions_from_csv(csv_file: str) -> List[Dict]:
    """
    Load accessions from a CSV file with columns 'accession' and optionally 'source'.
    """
    try:
        logger.info(f"Loading accessions from CSV file: {csv_file}")
        df = pd.read_csv(csv_file)
        
        if 'accession' not in df.columns:
            logger.error(f"CSV file must have 'accession' column")
            return []
        
        accessions = df['accession'].tolist()
        logger.info(f"Loaded {len(accessions)} accessions from CSV")
        return accessions
    except Exception as e:
        logger.error(f"Error loading accessions from CSV: {str(e)}")
        return []

def build_contaminant_database(output_fasta: str, accessions: List[str], 
                              max_sequences_per_batch: int = 10) -> List[SeqRecord]:
    """
    Build a comprehensive database of viral contaminants.
    Downloads sequences in batches for better reliability.
    """
    all_sequences = []
    
    # Process in batches to avoid overwhelming the NCBI API
    for i in range(0, len(accessions), max_sequences_per_batch):
        batch = accessions[i:i + max_sequences_per_batch]
        logger.info(f"Processing batch {i//max_sequences_per_batch + 1}/{(len(accessions)-1)//max_sequences_per_batch + 1} ({len(batch)} accessions)")
        
        for accession in batch:
            record = download_sequence_by_accession(accession)
            if record:
                # Store the accession for traceability
                record.annotations['original_accession'] = accession
                all_sequences.append(record)
            
        # Add a small delay between batches to be nice to NCBI servers
        if i + max_sequences_per_batch < len(accessions):
            time.sleep(1)
    
    # Write sequences to FASTA file
    if all_sequences:
        SeqIO.write(all_sequences, output_fasta, "fasta")
        logger.info(f"Wrote {len(all_sequences)} sequences to {output_fasta}")
    else:
        logger.warning("No sequences were successfully downloaded")
    
    return all_sequences

def create_metadata_file(sequences: List[SeqRecord], output_file: str) -> None:
    """Create a TSV file with metadata for all sequences."""
    metadata = []
    for record in sequences:
        entry = {
            "accession": record.id,
            "original_accession": record.annotations.get('original_accession', record.id),
            "description": record.description,
            "length": len(record.seq),
            "taxonomy": ";".join(record.annotations.get('taxonomy', [])),
            "organism": record.annotations.get('organism', ''),
            "source": record.annotations.get('source', '')
        }
        metadata.append(entry)
    
    df = pd.DataFrame(metadata)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Metadata written to {output_file}")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Build a viral contaminant database for decontamination pipelines")
    parser.add_argument("--accession-csv", required=True, 
                        help="CSV file with accessions (must have 'accession' column)")
    parser.add_argument("--output-dir", required=True, 
                        help="Output directory where files will be saved")
    parser.add_argument("--output-name", default="contaminant_db",
                        help="Base name for output files (default: contaminant_db)")
    parser.add_argument("--email", default="user@example.com", 
                        help="Email for NCBI Entrez queries")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Maximum number of sequences to download in each batch")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up file paths
    output_fasta = os.path.join(args.output_dir, f"{args.output_name}.fasta")
    output_metadata = os.path.join(args.output_dir, f"{args.output_name}_metadata.tsv")
    log_file = os.path.join(args.output_dir, f"{args.output_name}_build.log")
    
    # Set up logging
    logger = setup_logging(log_file)
    logger.info(f"Starting contaminant database build process")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Output FASTA: {output_fasta}")
    
    # Set up Entrez
    set_entrez_email(args.email)
    
    # Load accessions from CSV
    accessions = load_accessions_from_csv(args.accession_csv)
    if not accessions:
        logger.error("No accessions found in CSV file. Exiting.")
        sys.exit(1)
    
    # Build the database
    sequences = build_contaminant_database(output_fasta, accessions, args.batch_size)
    
    # Create metadata file
    if sequences:
        create_metadata_file(sequences, output_metadata)
        
    logger.info(f"Database build process completed")
    logger.info(f"Results saved to {args.output_dir}")
    logger.info(f"- FASTA database: {output_fasta}")
    logger.info(f"- Metadata file: {output_metadata}")
    logger.info(f"- Log file: {log_file}")


# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/build_contaminant_db.py \
#             --accession-csv /fs/project/PAS1117/ricardo/ssDNA_tool/DB/decont_accesion_list.csv \
#             --output-dir /fs/project/PAS1117/ricardo/ssDNA_tool/DB \    
#             --email pavan.4@osu.edu \
#             --batch-size 10