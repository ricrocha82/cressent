#!/usr/bin/env python3

import os
import argparse
import subprocess
import logging
from typing import List, Union, Optional
from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import sys

class CommandExecutor:
    @staticmethod
    def run_command(
        command: Union[str, List[str]],
        error_message: str,
        capture_output: bool = True,
        log_output: bool = True
    ) -> subprocess.CompletedProcess:
        """
        Enhanced helper function to run a shell command and handle errors.
        
        Args:
            command: Command to run (string or list)
            error_message: Error message to display if command fails
            capture_output: Whether to capture command output
            log_output: Whether to log command output
            
        Returns:
            CompletedProcess instance
        """
        try:
            if isinstance(command, list):
                cmd_str = " ".join(str(x) for x in command)
            else:
                cmd_str = command
                
            logging.info(f"Running command: {cmd_str}")
            
            result = subprocess.run(
                command,
                shell=isinstance(command, str),
                check=True,
                text=True,
                capture_output=capture_output
            )
            
            if log_output and result.stdout:
                logging.info(result.stdout)
                
            return result
            
        except subprocess.CalledProcessError as e:
            logging.error(f"{error_message}: {e}")
            if e.output:
                logging.error(f"Command output: {e.output}")
            if e.stderr:
                logging.error(f"Error output: {e.stderr}")
            raise RuntimeError(f"{error_message}: {e}")

class SequenceClusterer:
    def __init__(self, fasta_file, output_dir=None, threads=32, min_ani=95.0, min_tcov=85.0, min_qcov=0.0):
        """
        Initialize the SequenceClusterer with improved type hints and validation.
        
        Args:
            fasta_file: Path to input FASTA file
            output_dir: Directory for output files (default: current directory)
            threads: Number of threads for BLAST
            evalue: Minimum number of expected hits of similar quality (score) 
            min_ani: Minimum average identity for clustering
            min_tcov: Minimum target coverage
            min_qcov: Minimum query coverage
        """

        ### INPUTS ###
        self.fasta_file = Path(fasta_file)
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.threads = threads
        # self.evalue = float(evalue)
        self.min_ani = float(min_ani)
        self.min_tcov = float(min_tcov)
        self.min_qcov = float(min_qcov)

        # Preprocess input FASTA file
        self.cleaned_fasta = self.output_dir / f"cleaned_{self.fasta_file.name}"
        self.preprocess_fasta()

        # Automatically set the path to the modules directory
        self.modules_dir = Path(__file__).resolve().parent  # Directory where the script is located
        self.anicalc_script = self.modules_dir / "anicalc.py"
        self.aniclust_script = self.modules_dir / "aniclust.py"

        # Validate scripts exist
        if not self.anicalc_script.exists():
            raise FileNotFoundError(f"anicalc.py not found in {self.modules_dir}")
        if not self.aniclust_script.exists():
            raise FileNotFoundError(f"aniclust.py not found in {self.modules_dir}")

        # Create temp directory inside the output directory
        self.temp_dir = self.output_dir / "temp"
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        ### OUTPUTS ###
        self.db_path = self.temp_dir / "blast_db"
        self.blast_output = self.output_dir / "blast_results.tsv"
        self.ani_output = self.output_dir / "ani_results.tsv"
        self.cluster_output = self.output_dir / "clusters.tsv"
        self.cluster_fasta = self.output_dir / "cluster_sequences.fa"
        
        # Validate inputs and create output directory
        self.validate_inputs()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Set up logging
        self._setup_logging()
        
    def _setup_logging(self):
        """Configure logging for the clustering pipeline."""
        log_file = self.output_dir / "clustering.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
    
    def validate_inputs(self):
        """Validate input files and parameters."""
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"Input FASTA file not found: {self.fasta_file}")
        # if not self.utils_dir.exists():
        #     raise FileNotFoundError(f"Utils directory not found: {self.utils_dir}")
        if self.threads < 1:
            raise ValueError("Number of threads must be positive")
        if not (0 <= self.min_ani <= 100):
            raise ValueError("min_ani must be between 0 and 100")
        if not (0 <= self.min_tcov <= 100):
            raise ValueError("min_tcov must be between 0 and 100")
        if not (0 <= self.min_qcov <= 100):
            raise ValueError("min_qcov must be between 0 and 100")

    def detect_sequence_type(self):
        """Detect whether the input FASTA file contains nucleotide or protein sequences."""
        nucleotide_chars = {"A", "T", "C", "G", "N"}
        protein_chars = set("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids

        nucleotide_count = 0
        protein_count = 0
        total_chars = 0

        with open(self.fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue  # Skip headers
                line = line.strip().upper()  # Convert to uppercase for uniformity
                total_chars += len(line)
                nucleotide_count += sum(1 for c in line if c in nucleotide_chars)
                protein_count += sum(1 for c in line if c in protein_chars)

        # Determine sequence type
        if total_chars == 0:
            raise ValueError(f"Empty FASTA file: {self.fasta_file}")

        nucleotide_ratio = nucleotide_count / total_chars
        if nucleotide_ratio > 0.90:
            logging.info("Detected nucleotide sequences (blastn will be used).")
            return "nucleotide"
        else:
            logging.info("Detected protein sequences (blastp will be used).")
            return "protein"    

    def preprocess_fasta(self):
        """Cleans the FASTA file by replacing spaces with underscores and ensuring unique sequence IDs."""
        seen_ids = defaultdict(int)  # Track duplicate IDs
        with open(self.fasta_file, "r") as infile, open(self.cleaned_fasta, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    original_id = line.strip()[1:]  # Remove ">"
                    clean_id = original_id.replace(" ", "_")  # Replace spaces
                    seen_ids[clean_id] += 1

                    # Make ID unique if duplicate exists
                    if seen_ids[clean_id] > 1:
                        clean_id = f"{clean_id}_{seen_ids[clean_id]}"

                    outfile.write(f">{clean_id}\n")
                else:
                    outfile.write(line)
        # Use cleaned FASTA for processing
        self.fasta_file = self.cleaned_fasta
        

    
    def run_blast_makedb(self):
        """Creates a BLAST database, automatically detecting sequence type."""
        self.sequence_type = self.detect_sequence_type()
        db_type = "nucl" if self.sequence_type == "nucleotide" else "prot"

        cmd = f"makeblastdb -in {self.fasta_file} -dbtype {db_type} -out {self.db_path}"
        CommandExecutor.run_command(cmd, "Failed to create BLAST database")

    def run_blast_search(self):
        """Runs BLAST search using either blastn or blastp, depending on sequence type."""
        self.sequence_type = self.detect_sequence_type()
        
        if self.sequence_type == "nucleotide":
            blast_cmd = f"blastn -query {self.fasta_file} -db {self.db_path} -out {self.blast_output} -outfmt '6 std qlen slen' -max_target_seqs 10000 -num_threads {self.threads}"
        else:
            blast_cmd = f"blastp -query {self.fasta_file} -db {self.db_path} -out {self.blast_output} -outfmt '6 std qlen slen' -evalue 0.00001 -num_threads {self.threads}"

        CommandExecutor.run_command(blast_cmd, f"{self.sequence_type} BLAST search failed")

    def run_anicalc(self):
        """Calculates ANI from BLAST results using anicalc.py."""
        cmd = f"python {self.anicalc_script} -i {self.blast_output} -o {self.ani_output}"
        CommandExecutor.run_command(cmd, "ANI calculation failed")

    def run_aniclust(self):
        """Clusters sequences using aniclust.py."""
        cmd = f"python {self.aniclust_script} --fna {self.fasta_file} --ani {self.ani_output} --out {self.cluster_output} --min_ani {self.min_ani} --min_tcov {self.min_tcov} --min_qcov {self.min_qcov}"
        CommandExecutor.run_command(cmd, "Sequence clustering failed")
        # Add column names to the output file
        df = pd.read_csv(self.cluster_output, sep='\t', header=None)
        df.columns = ['Representative_Sequence', 'Sequences']
        contig_ids_to_keep = set(df['Representative_Sequence'].tolist())
        df.to_csv(self.cluster_output, sep='\t', index=False)

        return contig_ids_to_keep

    def create_representative_fasta_file(self, contig_ids_to_keep):
        print(f"\nCreating the FASTA file containing the representative sequences from the clustering step...")
        # Read all lines from the fasta file
        fasta_lines = []
        with open(self.fasta_file, 'r') as f:
            fasta_lines = f.readlines()

        # Filter the fasta file
        representative_fasta_file = self.cluster_fasta
        with open(representative_fasta_file, 'w') as f:
            write_line = False
            for line in fasta_lines:
                if line.startswith('>'):
                    contig_id = line.strip()[1:]
                    if contig_id in contig_ids_to_keep:
                        write_line = True
                        f.write(line)
                    else:
                        write_line = False
                else:
                    if write_line:
                        f.write(line)

        return representative_fasta_file


    def cleanup_files(self):
        """Remove all temporary files generated during the pipeline."""
        logging.info("Cleaning up temporary files...")

        # Delete all files in temp_dir
        for file in self.temp_dir.glob("*"):
            file.unlink()
            logging.info(f"Deleted: {file}")

        # Remove temp directory itself
        self.temp_dir.rmdir()
        logging.info(f"Deleted temp directory: {self.temp_dir}")


    def run_pipeline(self):
        """Executes the full clustering pipeline with proper cleanup."""
        try:
            logging.info("Starting sequence clustering pipeline")
            self.run_blast_makedb()
            self.run_blast_search()  # Automatically picks blastn or blastp
            self.run_anicalc()
            contig_ids_to_keep = self.run_aniclust()  # Capture the returned set
            self.create_representative_fasta_file(contig_ids_to_keep)
            logging.info(f"Clustering completed successfully! Results saved in {self.cluster_output}")
        except Exception as e:
            logging.error(f"Pipeline failed: {e}")
            raise
        finally:
            self.cleanup_files()  # Cleanup temporary files after execution

def main():
    parser = argparse.ArgumentParser(
        description="Sequence clustering using BLAST, anicalc, and aniclust.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input_fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", default=".", help="Directory to save output")
    parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for BLAST")
    parser.add_argument("--min_ani", type=float, default=95.0, help="Minimum average identity for clustering")
    parser.add_argument("--min_tcov", type=float, default=85.0, help="Minimum target coverage")
    parser.add_argument("--min_qcov", type=float, default=0.0, help="Minimum query coverage")
    
    args = parser.parse_args()

    output_dir = args.output

    # Determine the input FASTA full path
    def validate_fasta(filename):
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                print("FASTA checked.")
                input_fasta = os.path.join(output_dir, filename)
                return input_fasta
            else:
                sys.exit("Error: Input file is not in the FASTA format.\n")
                logging.info(f"Using: {input_fasta} is not in the FASTA format")
    # check fasta
    input_fasta = validate_fasta(args.input_fasta)
    
    try:
        clusterer = SequenceClusterer(
            fasta_file=input_fasta,
            output_dir=output_dir,
            threads=args.threads,
            # evalue=args.evalue,
            min_ani=args.min_ani,
            min_tcov=args.min_tcov,
            min_qcov=args.min_qcov
        )
        clusterer.run_pipeline()
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        exit(1)

if __name__ == "__main__":
    main()


## run with blastn
# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/cluster.py \
#           -i /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/1.virus_detection/vs2_1.5kb/ssDNA/all_ssDNA.fasta \
#           -o /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/output_clusters

# output /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/output_clusters/clusters.tsv is the file with the first column as the cluster representative

## run with blastp 
# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/cluster.py \
#           -i /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/reps.fa \
#           -o /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/output_clusters
