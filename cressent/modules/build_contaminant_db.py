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

def download_sequence_and_proteins(accession: str, retries: int = 3) -> tuple:
    """
    Download a sequence from NCBI using its accession number and extract its protein sequences.
    Includes retry logic for reliability.
    
    Returns:
        Tuple of (nucleotide_record, list_of_protein_records)
    """
    for attempt in range(retries):
        try:
            logger.info(f"Downloading sequence: {accession}")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Validate that the sequence has content
            if not record.seq or len(record.seq) == 0 or str(record.seq) == '':
                logger.warning(f"Downloaded empty sequence for {accession}")
                if attempt < retries - 1:
                    logger.info(f"Retrying {accession}...")
                    time.sleep(2)
                    continue
                else:
                    return None, []
            
            # Add the accession to the description for better tracking
            if not record.description:
                record.description = f"Viral contaminant {accession}"
            else:
                record.description = f"{record.description} [Contaminant]"
            
            # Extract protein IDs from CDS features
            protein_ids = []
            for feature in record.features:
                if feature.type == "CDS":
                    if "protein_id" in feature.qualifiers:
                        protein_ids.append(feature.qualifiers["protein_id"][0])
            
            # Download protein sequences
            proteins = []
            if protein_ids:
                logger.info(f"Found {len(protein_ids)} protein IDs for {accession}")
                for protein_id in protein_ids:
                    try:
                        logger.info(f"Downloading protein: {protein_id} from {accession}")
                        protein_handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
                        protein_record = SeqIO.read(protein_handle, "genbank")
                        protein_handle.close()
                        
                        # Add the original accession to annotations
                        protein_record.annotations['original_accession'] = accession
                        # Modify the description to indicate this is from a contaminant
                        if not protein_record.description:
                            protein_record.description = f"Protein from viral contaminant {accession}"
                        else:
                            protein_record.description = f"{protein_record.description} [Contaminant]"
                        
                        proteins.append(protein_record)
                    except Exception as e:
                        logger.warning(f"Failed to download protein {protein_id}: {str(e)}")
                        # Continue with other proteins even if one fails
                        continue
            else:
                logger.info(f"No protein IDs found for {accession}")
                
            return record, proteins
            
        except Exception as e:
            logger.warning(f"Attempt {attempt+1}/{retries} failed for {accession}: {str(e)}")
            if attempt < retries - 1:
                time.sleep(2)  # Wait 2 seconds before retrying
            else:
                logger.error(f"Failed to download {accession} after {retries} attempts")
                return None, []

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
                              max_sequences_per_batch: int = 10,
                              output_protein_fasta: str = None) -> tuple:
    """
    Build a comprehensive database of viral contaminants.
    Downloads sequences in batches for better reliability.
    
    Returns:
        Tuple of (nucleotide_sequences, protein_sequences)
    """
    all_sequences = []
    all_proteins = []
    
    # First make sure accessions is actually a list
    if not isinstance(accessions, list):
        logger.error(f"Accessions must be a list, got {type(accessions)}")
        return [], []
    
    # Process in batches to avoid overwhelming the NCBI API
    for i in range(0, len(accessions), max_sequences_per_batch):
        batch = accessions[i:i + max_sequences_per_batch]
        logger.info(f"Processing batch {i//max_sequences_per_batch + 1}/{(len(accessions)-1)//max_sequences_per_batch + 1} ({len(batch)} accessions)")
        
        for accession in batch:
            # Download nucleotide sequence and protein sequences
            record, proteins = download_sequence_and_proteins(accession)
            
            # Process nucleotide sequence
            if record:
                # Store the accession for traceability
                record.annotations['original_accession'] = accession
                
                # Verify that the sequence has content before adding it
                if record.seq and len(record.seq) > 0 and not hasattr(record.seq, 'undefined') and str(record.seq) != '':
                    all_sequences.append(record)
                else:
                    logger.warning(f"Sequence for {accession} has no content and will be skipped")
            
            # Process protein sequences
            if proteins:
                # Verify each protein sequence before adding
                for protein in proteins:
                    if protein.seq and len(protein.seq) > 0 and not hasattr(protein.seq, 'undefined') and str(protein.seq) != '':
                        all_proteins.append(protein)
                    else:
                        logger.warning(f"Protein sequence for {protein.id} has no content and will be skipped")
            
        # Add a small delay between batches to be nice to NCBI servers
        if i + max_sequences_per_batch < len(accessions):
            time.sleep(1)
    
    # Write nucleotide sequences to FASTA file
    if all_sequences:
        # Additional check before writing to file
        valid_sequences = []
        for seq in all_sequences:
            try:
                # Test if we can convert the sequence to string
                str(seq.seq)
                valid_sequences.append(seq)
            except Exception as e:
                logger.warning(f"Skipping sequence {seq.id}: {str(e)}")
        
        SeqIO.write(valid_sequences, output_fasta, "fasta")
        logger.info(f"Wrote {len(valid_sequences)} nucleotide sequences to {output_fasta}")
    else:
        logger.warning("No nucleotide sequences were successfully downloaded")
    
    # Write protein sequences to FASTA file if output file is specified
    if all_proteins and output_protein_fasta:
        # Additional check before writing to file
        valid_proteins = []
        for protein in all_proteins:
            try:
                # Test if we can convert the sequence to string
                str(protein.seq)
                valid_proteins.append(protein)
            except Exception as e:
                logger.warning(f"Skipping protein {protein.id}: {str(e)}")
        
        SeqIO.write(valid_proteins, output_protein_fasta, "fasta")
        logger.info(f"Wrote {len(valid_proteins)} protein sequences to {output_protein_fasta}")
    else:
        if output_protein_fasta:
            logger.warning("No protein sequences were successfully downloaded")
    
    return (all_sequences, all_proteins)

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
    parser.add_argument("-o","--output", required=True, 
                        help="Path to the output directory (Default: working directory)")
    parser.add_argument("--output-name", default="contaminant_db",
                        help="Base name for output files (default: contaminant_db)")
    parser.add_argument("--email", default="user@example.com", 
                        help="Email for NCBI Entrez queries (not required, default: user@example.com)")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Maximum number of sequences to download in each batch (default = 10)")
    return parser.parse_args()

def main():
    """Main function to run the contaminant database builder (for CLI import)."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Set up file paths
    output_fasta = os.path.join(args.output, f"{args.output_name}.fasta")
    output_protein_fasta = os.path.join(args.output, f"{args.output_name}_proteins.fasta")
    output_metadata = os.path.join(args.output, f"{args.output_name}_metadata.tsv")
    output_protein_metadata = os.path.join(args.output, f"{args.output_name}_protein_metadata.tsv")
    log_file = os.path.join(args.output, f"{args.output_name}_build.log")
    
    # Set up global logger
    global logger
    logger = setup_logging(log_file)
    logger.info(f"Starting contaminant database build process")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Output FASTA: {output_fasta}")
    logger.info(f"Output protein FASTA: {output_protein_fasta}")
    
    # Set up Entrez
    set_entrez_email(args.email)
    
    # Load accessions from CSV
    accessions = load_accessions_from_csv(args.accession_csv)
    if not accessions:
        logger.error("No accessions found in CSV file. Exiting.")
        return 1
        
    logger.info(f"Loaded {len(accessions)} accessions from CSV")
    
    # Build the database
    sequences, proteins = build_contaminant_database(
        output_fasta=output_fasta,
        accessions=accessions,
        max_sequences_per_batch=args.batch_size,
        output_protein_fasta=output_protein_fasta
    )
    
    # Create metadata files
    if sequences:
        create_metadata_file(sequences, output_metadata)
    
    if proteins:
        create_metadata_file(proteins, output_protein_metadata)
        
    logger.info(f"Database build process completed")
    logger.info(f"Results saved to {args.output}")
    logger.info(f"- FASTA database: {output_fasta}")
    logger.info(f"- Protein FASTA database: {output_protein_fasta}")
    logger.info(f"- Metadata file: {output_metadata}")
    logger.info(f"- Protein metadata file: {output_protein_metadata}")
    logger.info(f"- Log file: {log_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

