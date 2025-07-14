#!/usr/bin/env python3
"""
Taxonomy-based Database Builder for ssDNA Tool
Builds custom databases using taxonomy selection from ICTV data
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from Bio import Entrez, SeqIO
from typing import List, Dict, Set
import time
import logging

# Configure logging
def setup_logging(output_dir: str, taxonomy_name: str) -> None:
    """Setup logging configuration with log file in output directory"""
    log_dir = os.path.join(output_dir, taxonomy_name)
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    log_file = os.path.join(log_dir, 'db_builder.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logging.info(f"Log file created at: {log_file}")
    
# After you have taxonomy_name and args.output_dir
def setup_directories(base_dir: str, taxonomy_name: str) -> Dict[str, str]:
    """Create necessary directories for the pipeline"""
    dirs = {
        'base': base_dir,
        'taxonomy': os.path.join(base_dir, taxonomy_name),
        'raw_aa': os.path.join(base_dir, taxonomy_name, 'raw_aa'),
        'cd_hit': os.path.join(base_dir, taxonomy_name, 'cd_hit'),
        'diamond': os.path.join(base_dir, taxonomy_name, 'diamond'),
        'mcl': os.path.join(base_dir, taxonomy_name, 'mcl'),
        'annotated': os.path.join(base_dir, taxonomy_name, 'annotated'),
        'unannotated': os.path.join(base_dir, taxonomy_name, 'unannotated')
    }
    
    for dir_path in dirs.values():
        Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    return dirs

def load_taxonomy_data(taxonomy_file: str) -> pd.DataFrame:
    """Load taxonomy data from CSV file"""
    try:
        df = pd.read_csv(taxonomy_file)
        logging.info(f"Loaded taxonomy data with {len(df)} records")
        return df
    except Exception as e:
        logging.error(f"Error loading taxonomy file: {e}")
        sys.exit(1)

def get_available_taxonomies(df: pd.DataFrame, taxonomy_level: str) -> List[str]:
    """Get unique taxonomies available at specified level"""
    if taxonomy_level not in df.columns:
        logging.error(f"Taxonomy level '{taxonomy_level}' not found in data")
        sys.exit(1)
    
    taxonomies = df[taxonomy_level].dropna().unique().tolist()
    logging.info(f"Found {len(taxonomies)} unique {taxonomy_level} taxonomies")
    return sorted(taxonomies)

def filter_by_taxonomy(df: pd.DataFrame, taxonomy_level: str, selected_taxonomies: List[str]) -> pd.DataFrame:
    """Filter dataframe by selected taxonomies"""
    filtered_df = df[df[taxonomy_level].isin(selected_taxonomies)]
    logging.info(f"Filtered to {len(filtered_df)} records for selected {taxonomy_level}")
    return filtered_df

def download_sequences_from_accessions(accessions: List[str], output_file: str, email: str, batch_size: int = 50) -> None:
    """Download protein sequences from NCBI using accession numbers"""
    Entrez.email = email
    
    # Remove duplicates and clean accessions
    unique_accessions = list(set([acc.strip() for acc in accessions if acc.strip()]))
    logging.info(f"Downloading sequences for {len(unique_accessions)} unique accessions")
    
    with open(output_file, 'w') as out_file:
        for i in range(0, len(unique_accessions), batch_size):
            batch = unique_accessions[i:i + batch_size]
            logging.info(f"Processing batch {i//batch_size + 1}/{(len(unique_accessions) + batch_size - 1)//batch_size}")
            
            try:
                # Search for nucleotide records
                for accession in batch:
                    try:
                        # First try to get nucleotide record
                        search_handle = Entrez.esearch(db="nuccore", term=accession, retmax=1)
                        search_results = Entrez.read(search_handle)
                        search_handle.close()
                        
                        if search_results["IdList"]:
                            nucleotide_id = search_results["IdList"][0]
                            
                            # Get feature table to extract protein IDs
                            fetch_handle = Entrez.efetch(db="nuccore", id=nucleotide_id, rettype="ft", retmode="text")
                            feature_table = fetch_handle.read()
                            fetch_handle.close()
                            
                            # Extract protein IDs
                            protein_ids = []
                            for line in feature_table.splitlines():
                                if "protein_id" in line:
                                    parts = line.split()
                                    if len(parts) > 1:
                                        protein_id = parts[1].strip('"')
                                        protein_ids.append(protein_id)
                            
                            # Download protein sequences
                            for protein_id in protein_ids:
                                try:
                                    protein_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                                    fasta_data = protein_handle.read()
                                    protein_handle.close()
                                    out_file.write(fasta_data)
                                    time.sleep(0.1)  # Rate limiting
                                except Exception as e:
                                    logging.warning(f"Failed to download protein {protein_id}: {e}")
                        else:
                            logging.warning(f"No nucleotide record found for {accession}")
                    
                    except Exception as e:
                        logging.warning(f"Error processing accession {accession}: {e}")
                    
                    time.sleep(0.1)  # Rate limiting
                
            except Exception as e:
                logging.error(f"Error in batch processing: {e}")
                continue

def run_cd_hit(input_file: str, output_file: str, identity: float = 0.95, coverage: float = 0.9, threads: int = 8) -> bool:
    """Run CD-HIT clustering"""

    Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)

    cmd = [
        'cd-hit',
        '-i', input_file,
        '-o', output_file,
        '-c', str(identity),
        '-aS', str(coverage),
        '-M', '16000',
        '-T', str(threads)
    ]
    
    try:
        logging.info(f"Running CD-HIT: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info("CD-HIT completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"CD-HIT failed: {e}")
        return False

def run_diamond(input_file: str, output_dir: str, threads: int = 24) -> bool:
    """Run Diamond BLASTP"""
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    db_file = os.path.join(output_dir, base_name)
    output_file = os.path.join(output_dir, f"{base_name}.txt")
    
    # Create database
    makedb_cmd = [
        'diamond', 'makedb',
        '--in', input_file,
        '-d', db_file,
        '--threads', str(threads)
    ]
    
    try:
        logging.info(f"Creating Diamond database: {' '.join(makedb_cmd)}")
        subprocess.run(makedb_cmd, check=True)
        
        # Run BLASTP
        blastp_cmd = [
            'diamond', 'blastp',
            '--db', db_file,
            '--query', input_file,
            '--out', output_file,
            '--outfmt', '6',
            '--threads', str(threads),
            '--evalue', '0.00001',
            '--masking', '0',
            '--sensitive'
        ]
        
        logging.info(f"Running Diamond BLASTP: {' '.join(blastp_cmd)}")
        subprocess.run(blastp_cmd, check=True)
        logging.info("Diamond BLASTP completed successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Diamond failed: {e}")
        return False

def run_mcl(diamond_output: str, mcl_output: str, inflation: float = 1.5, threads: int = 24) -> bool:
    """Run MCL clustering"""
    # Prepare MCL input (take columns 1, 2, 11 from diamond output)
    mcl_input = diamond_output.replace('.txt', '_mcl_input')
    
    try:
        # Prepare input file
        with open(diamond_output, 'r') as infile, open(mcl_input, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 11:
                    outfile.write(f"{parts[0]}\t{parts[1]}\t{parts[10]}\n")
        
        # Run MCL
        mcl_cmd = [
            'mcl', mcl_input,
            '--abc',
            '-I', str(inflation),
            '-te', str(threads),
            '-o', mcl_output,
            '--abc-neg-log10'
        ]
        
        logging.info(f"Running MCL: {' '.join(mcl_cmd)}")
        subprocess.run(mcl_cmd, check=True)
        logging.info("MCL completed successfully")
        
        # Clean up temporary file
        os.remove(mcl_input)
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"MCL failed: {e}")
        return False

def process_mcl_output(mcl_file: str, fasta_file: str, output_file: str) -> None:
    """Process MCL output and add protein descriptions"""
    # Read MCL clusters
    clusters = []
    with open(mcl_file, 'r') as f:
        for i, line in enumerate(f):
            proteins = line.strip().split('\t')
            for protein in proteins:
                if protein:
                    clusters.append({'cluster': f'cl_{i+1}', 'protein_id': protein})
    
    df_clusters = pd.DataFrame(clusters)
    
    # Read FASTA file to get descriptions
    fasta_data = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        fasta_data.append({
            'protein_id': record.id,
            'protein_description': record.description,
            'protein_sequence': str(record.seq)
        })
    
    df_fasta = pd.DataFrame(fasta_data)
    
    # Merge clusters with FASTA data
    df_merged = df_clusters.merge(df_fasta, on='protein_id', how='left')
    
    # Save results
    df_merged.to_csv(output_file, index=False)
    logging.info(f"Processed MCL output saved to {output_file}")

def annotate_proteins(df: pd.DataFrame, output_dirs: Dict[str, str]) -> None:
    """Annotate proteins into caps, reps, and ORFs"""
    
    # Capsid proteins
    df_caps = df[
        df['protein_description'].str.contains('cap|Cap|capsid|Capsid|VP|coat', case=False, na=False) &
        ~df['protein_description'].str.contains('non-capsid|non-Capsid', case=False, na=False)
    ]
    
    if not df_caps.empty:
        # Get clusters with at least 8 capsid proteins
        cluster_counts = df_caps['cluster'].value_counts()
        clusters_to_keep = cluster_counts[cluster_counts >= 8].index.tolist()
        
        if not clusters_to_keep:
            clusters_to_keep = cluster_counts[cluster_counts >= 1].index.tolist()
        
        if clusters_to_keep:
            df_caps_filtered = df_caps[df_caps['cluster'].isin(clusters_to_keep)]
            
            # Create caps directory
            caps_dir = os.path.join(output_dirs['annotated'], 'caps')
            Path(caps_dir).mkdir(exist_ok=True)
            
            # Save CSV
            df_caps_filtered[['cluster', 'protein_id', 'protein_description']].to_csv(
                os.path.join(caps_dir, 'caps.csv'), index=False
            )
            
            # Save FASTA files by cluster
            for cluster in clusters_to_keep:
                cluster_seqs = df_caps_filtered[df_caps_filtered['cluster'] == cluster]
                fasta_file = os.path.join(caps_dir, f'{cluster}.fa')
                with open(fasta_file, 'w') as f:
                    for _, row in cluster_seqs.iterrows():
                        f.write(f">{row['protein_id']} {row['protein_description']}\n")
                        f.write(f"{row['protein_sequence']}\n")
    
    # Rep proteins
    df_reps = df[
        df['protein_description'].str.contains('rep|Rep|replication|Replication|replicase|Replicase', case=False, na=False) &
        ~df['protein_description'].str.contains('non-rep|non-Rep|non-replication|non-Replication|non-replicase|non-Replicase', case=False, na=False)
    ]
    
    if not df_reps.empty:
        cluster_counts = df_reps['cluster'].value_counts()
        clusters_to_keep = cluster_counts[cluster_counts >= 1].index.tolist()
        
        if clusters_to_keep:
            df_reps_filtered = df_reps[df_reps['cluster'].isin(clusters_to_keep)]
            
            # Create reps directory
            reps_dir = os.path.join(output_dirs['annotated'], 'reps')
            Path(reps_dir).mkdir(exist_ok=True)
            
            # Save CSV
            df_reps_filtered[['cluster', 'protein_id', 'protein_description']].to_csv(
                os.path.join(reps_dir, 'reps.csv'), index=False
            )
            
            # Save FASTA files by cluster
            for cluster in clusters_to_keep:
                cluster_seqs = df_reps_filtered[df_reps_filtered['cluster'] == cluster]
                fasta_file = os.path.join(reps_dir, f'{cluster}.fa')
                with open(fasta_file, 'w') as f:
                    for _, row in cluster_seqs.iterrows():
                        f.write(f">{row['protein_id']} {row['protein_description']}\n")
                        f.write(f"{row['protein_sequence']}\n")
    
    # Get clusters used for caps and reps
    used_clusters = set()
    if not df_caps.empty:
        used_clusters.update(df_caps['cluster'].unique())
    if not df_reps.empty:
        used_clusters.update(df_reps['cluster'].unique())
    
    # ORFs (unannotated)
    df_orfs = df[~df['cluster'].isin(used_clusters)]
    
    if not df_orfs.empty:
        cluster_counts = df_orfs['cluster'].value_counts()
        clusters_to_keep = cluster_counts[cluster_counts >= 10].index.tolist()
        
        if clusters_to_keep:
            df_orfs_filtered = df_orfs[df_orfs['cluster'].isin(clusters_to_keep)]
            
            # Save CSV
            df_orfs_filtered[['cluster', 'protein_id', 'protein_description']].to_csv(
                os.path.join(output_dirs['unannotated'], 'orfs.csv'), index=False
            )
            
            # Save FASTA files by cluster
            for cluster in clusters_to_keep:
                cluster_seqs = df_orfs_filtered[df_orfs_filtered['cluster'] == cluster]
                fasta_file = os.path.join(output_dirs['unannotated'], f'{cluster}.fa')
                with open(fasta_file, 'w') as f:
                    for _, row in cluster_seqs.iterrows():
                        f.write(f">{row['protein_id']} {row['protein_description']}\n")
                        f.write(f"{row['protein_sequence']}\n")

def main():
    parser = argparse.ArgumentParser(description='Build taxonomy-based database for ssDNA tool')
    parser.add_argument('-t', '--taxonomy-file', required=True, 
                       help='Path to taxonomy_accession_number.csv file')
    parser.add_argument('-l', '--taxonomy-level', required=True,
                       choices=['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 
                               'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 
                               'Genus', 'Subgenus', 'Species'],
                       help='Taxonomy level to use for selection')
    parser.add_argument('-s', '--selected-taxonomies', nargs='+', 
                       help='Selected taxonomies (if not provided, will list available options)')
    parser.add_argument('-o', '--output-dir', default=".",
                       help='Output directory for database')
    parser.add_argument('-e', '--email', default="my_email@gmail.com",
                       help='Email for NCBI Entrez')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads to use (default: 8)')
    parser.add_argument('--cd-hit-identity', type=float, default=0.95,
                       help='CD-HIT identity threshold (default: 0.95)')
    parser.add_argument('--mcl-inflation', type=float, default=1.5,
                       help='MCL inflation parameter (default: 1.5)')
    
    args = parser.parse_args()

    # Create taxonomy name for output
    taxonomy_name = "_".join(args.selected_taxonomies).replace(" ", "_")

    # After you have taxonomy_name and args.output_dir
    setup_logging(args.output_dir, taxonomy_name)
    
    # Load taxonomy data
    df_taxonomy = load_taxonomy_data(args.taxonomy_file)
    
    # If no taxonomies selected, show available options
    if not args.selected_taxonomies:
        available_taxonomies = get_available_taxonomies(df_taxonomy, args.taxonomy_level)
        print(f"\nAvailable {args.taxonomy_level} taxonomies:")
        for i, tax in enumerate(available_taxonomies, 1):
            print(f"{i:3d}. {tax}")
        print(f"\nTotal: {len(available_taxonomies)} taxonomies")
        print(f"\nUsage: python {sys.argv[0]} -t {args.taxonomy_file} -l {args.taxonomy_level} -s TAXONOMY1 TAXONOMY2 ... -o {args.output_dir} -e {args.email}")
        return
    
    # Filter by selected taxonomies
    df_filtered = filter_by_taxonomy(df_taxonomy, args.taxonomy_level, args.selected_taxonomies)
    
    if df_filtered.empty:
        logging.error("No records found for selected taxonomies")
        return
    
    # Setup directories
    dirs = setup_directories(args.output_dir, taxonomy_name)
    
    # Get accession numbers
    accessions = df_filtered['Virus GENBANK accession'].dropna().tolist()
    
    if not accessions:
        logging.error("No valid accession numbers found")
        return
    
    # Step 1: Download sequences
    raw_fasta = os.path.join(dirs['raw_aa'], f"{taxonomy_name}.fa")
    logging.info("Step 1: Downloading sequences from NCBI")
    download_sequences_from_accessions(accessions, raw_fasta, args.email)
    
    # Step 2: CD-HIT clustering
    clustered_fasta = os.path.join(dirs['cd_hit'], f"{taxonomy_name}.fa")
    logging.info("Step 2: Running CD-HIT clustering")
    if not run_cd_hit(raw_fasta, clustered_fasta, args.cd_hit_identity, threads=args.threads):
        logging.error("CD-HIT failed")
        return
    
    # Step 3: Diamond BLASTP
    logging.info("Step 3: Running Diamond BLASTP")
    if not run_diamond(clustered_fasta, dirs['diamond'], threads=args.threads):
        logging.error("Diamond failed")
        return
    
    # Step 4: MCL clustering
    diamond_output = os.path.join(dirs['diamond'], f"{taxonomy_name}.txt")
    mcl_output = os.path.join(dirs['mcl'], f"{taxonomy_name}_mcl_evalue")
    logging.info("Step 4: Running MCL clustering")
    if not run_mcl(diamond_output, mcl_output, args.mcl_inflation, threads=args.threads):
        logging.error("MCL failed")
        return
    
    # Step 5: Process MCL output
    protein_description_file = os.path.join(dirs['mcl'], f"{taxonomy_name}_protein_description.csv")
    logging.info("Step 5: Processing MCL output")
    process_mcl_output(mcl_output, clustered_fasta, protein_description_file)
    
    # Step 6: Annotate proteins
    logging.info("Step 6: Annotating proteins")
    df_proteins = pd.read_csv(protein_description_file)
    annotate_proteins(df_proteins, dirs)
    
    logging.info(f"Database building completed successfully!")
    logging.info(f"Results saved in: {dirs['taxonomy']}")

if __name__ == "__main__":
    main()


