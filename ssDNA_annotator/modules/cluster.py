#!/usr/bin/env python3

import os
import argparse
import subprocess
import logging
from typing import List, Union, Optional
from pathlib import Path

UTILS_DIR = "/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/utils"

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
    def __init__(
        self,
        fasta_file: Union[str, Path],
        output_dir: Optional[Union[str, Path]] = None,
        utils_dir: Optional[Union[str, Path]] = UTILS_DIR,
        threads: int = 32,
        # evalue: float = 0.00001,
        min_ani: float = 95.0,
        min_tcov: float = 85.0,
        min_qcov: float = 0.0
    ):
        """
        Initialize the SequenceClusterer with improved type hints and validation.
        
        Args:
            fasta_file: Path to input FASTA file
            output_dir: Directory for output files (default: current directory)
            utils_dir: Directory containing utility scripts
            threads: Number of threads for BLAST
            evalue: Minimum number of expected hits of similar quality (score) 
            min_ani: Minimum ANI for clustering
            min_tcov: Minimum target coverage
            min_qcov: Minimum query coverage
        """
        self.fasta_file = Path(fasta_file)
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.utils_dir = Path(utils_dir)
        self.threads = threads
        # self.evalue = float(evalue)
        self.min_ani = float(min_ani)
        self.min_tcov = float(min_tcov)
        self.min_qcov = float(min_qcov)
        
        # Define output file paths
        self.db_path = self.output_dir / "blast_db"
        self.blast_output = self.output_dir / "blast_results.tsv"
        self.ani_output = self.output_dir / "ani_results.tsv"
        self.cluster_output = self.output_dir / "clusters.tsv"
        
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
        if not self.utils_dir.exists():
            raise FileNotFoundError(f"Utils directory not found: {self.utils_dir}")
        if self.threads < 1:
            raise ValueError("Number of threads must be positive")
        if not (0 <= self.min_ani <= 100):
            raise ValueError("min_ani must be between 0 and 100")
        if not (0 <= self.min_tcov <= 100):
            raise ValueError("min_tcov must be between 0 and 100")
        if not (0 <= self.min_qcov <= 100):
            raise ValueError("min_qcov must be between 0 and 100")

    def run_blast_makedb(self):
        """Creates a BLAST database from the input FASTA file."""

        cmd = f"makeblastdb -in {self.fasta_file} -dbtype nucl -out {self.db_path}"
        CommandExecutor.run_command(cmd, "Failed to create BLAST database")

    def run_blast_blastp(self):
        """Runs blastn for all-vs-all similarity search with optimized parameters."""

        cmd = f"blastn -query {self.fasta_file} -db {self.db_path} -out {self.blast_output} -outfmt '6 std qlen slen' -max_target_seqs 10000 -num_threads {self.threads}"
        CommandExecutor.run_command(cmd, "blastn search failed")

    def run_anicalc(self):
        """Calculates ANI from BLAST results using anicalc.py."""
        cmd = f"python {self.utils_dir / "anicalc.py"} -i {self.blast_output} -o {self.ani_output}"
        CommandExecutor.run_command(cmd, "ANI calculation failed")

    def run_aniclust(self):
        """Clusters sequences using aniclust.py."""
        cmd = f"python {self.utils_dir / "aniclust.py"} --fna {self.fasta_file} --ani {self.ani_output} --out {self.cluster_output} --min_ani {self.min_ani} --min_tcov {self.min_tcov} --min_qcov {self.min_qcov}"
        CommandExecutor.run_command(cmd, "Sequence clustering failed")

    def run_pipeline(self):
        """Executes the full clustering pipeline with proper cleanup."""
        try:
            logging.info("Starting sequence clustering pipeline")
            self.run_blast_makedb()
            self.run_blast_blastp()
            self.run_anicalc()
            self.run_aniclust()
            logging.info(f"Clustering completed successfully! Results saved in {self.cluster_output}")
        except Exception as e:
            logging.error(f"Pipeline failed: {e}")
            raise
        finally:
            # Cleanup temporary files
            if self.db_path.exists():
                self.db_path.unlink()

def main():
    parser = argparse.ArgumentParser(
        description="Sequence clustering using BLAST, anicalc, and aniclust.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", default=".", help="Directory to save output")
    parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for BLAST")
    parser.add_argument("--min_ani", type=float, default=95.0, help="Minimum ANI for clustering")
    parser.add_argument("--min_tcov", type=float, default=85.0, help="Minimum target coverage")
    parser.add_argument("--min_qcov", type=float, default=0.0, help="Minimum query coverage")
    
    args = parser.parse_args()
    
    try:
        clusterer = SequenceClusterer(
            fasta_file=args.input,
            output_dir=args.output,
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


## ONLY Nucleotide 
# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/cluster.py \
#           -i /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/1.virus_detection/vs2_1.5kb/ssDNA/all_ssDNA.fasta \
#           -o /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/output_clusters

# output /fs/project/PAS1117/ricardo/ssDNA_test/test_data/virus/output_clusters/clusters.tsv is the file with the first column as the cluster representative