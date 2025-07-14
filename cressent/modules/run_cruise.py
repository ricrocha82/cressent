#!/usr/bin/env python3

# this is a modified version of Cruise (CRUcivirus Iteron SEarch) searches cressDNA virus genomes for iterons from 
# Jones, A., Kasun, G.W., Stover, J., Stedman, K.M. and de la Higuera, I., 2023. CRUISE, a Tool for the Detection of Iterons in Circular Rep-Encoding Single-Stranded DNA Viruses. Microbiology Resource Announcements, 12(1), pp.e01123-22.
# https://doi.org/10.1128/mra.01123-22
# CIte if using this module
# All rights reserved.

import os
import time
import logging
from pathlib import Path
import shutil

# Import from local module files
from . import cruise
from . import args

from Bio import SeqIO
import sys

def setup_logger(output_dir):
    """Set up logging configuration with file output in the specified directory"""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    # Create formatters and handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler
    log_file = Path(output_dir) / 'cruise.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    return logger

def input_convert(input_fasta_file, input_gff, input_folder):
    """Convert input files to individual GFF files"""
    fasta_dict = {}
    
    # Read FASTA file
    with open(input_fasta_file) as input_fasta:
        current_name = None
        current_sequence = []
        
        for line in input_fasta:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    fasta_dict[current_name] = "".join(current_sequence)
                current_name = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)
        
        if current_name:
            fasta_dict[current_name] = "".join(current_sequence)

    # Process GFF file and create individual files
    input_folder = Path(input_folder)
    input_folder.mkdir(parents=True, exist_ok=True)
    
    with open(input_gff) as gff:
        for name, sequence in fasta_dict.items():
            output_file = input_folder / f"{name}.gff"
            with output_file.open('w') as fout:
                # Write GFF3 header
                fout.write(f"""##gff-version 3
                            ##source-version geneious 2020.1.2
                            ##sequence-region\t{name}\t1\t{len(sequence)+1}
                            {name}\tGeneious\tregion\t1\t{len(sequence)+1}\t.\t+\t0\tIs_circular=true
                            """)
                
                # Write matching GFF entries
                gff.seek(0)
                for line in gff:
                    if name in line:
                        fout.write(line.replace("   ", "\t"))
                
                # Write sequence
                fout.write(f"##FASTA\n>{name}\n{sequence}\n")

def output_convert(output_folder, output_gff, output_annotations):
    """Combine output files into single GFF"""
    annotation_sources = output_annotations.split()
    position_index = int(annotation_sources[0])
    annotation_sources = annotation_sources[1:]
    
    # Validate position index
    filter_annotations = True
    if position_index == -1 or position_index not in range(1, 10):
        filter_annotations = False
        logging.warning("Invalid position index, preserving all annotations")

    output_folder = Path(output_folder)
    with open(output_gff, 'w') as output:
        for entry in output_folder.glob('*.gff'):
            with entry.open() as fin:
                in_sequence = False
                for line in fin:
                    if "##FASTA" in line:
                        break
                    elif "##sequence-region" in line:
                        in_sequence = True
                    elif in_sequence:
                        if filter_annotations:
                            fields = line.split()
                            feature = fields[position_index-1]
                            if feature in annotation_sources:
                                output.write(line)
                        else:
                            output.write(line)

def cleanup_internal_dirs(internal_dir):
    """
    Clean up internal directories recursively.
    Args:
        internal_dir (Path or str): Path to the internal directory
    """
    internal_dir = Path(internal_dir)
    
    try:
        if internal_dir.exists():
            # Use shutil.rmtree to recursively remove directory and contents
            shutil.rmtree(internal_dir)
            logging.info(f"Cleaned up internal working directory: {internal_dir}")
            
    except PermissionError as e:
        logging.error(f"Permission denied when trying to remove directory {internal_dir}: {e}")
    except FileNotFoundError as e:
        logging.warning(f"Directory {internal_dir} already removed or not found: {e}")
    except Exception as e:
        logging.error(f"Unexpected error when cleaning up directory {internal_dir}: {e}")
        # Try alternative cleanup if shutil.rmtree fails
        try:
            # Remove files first
            for file_path in internal_dir.rglob('*'):
                if file_path.is_file():
                    file_path.unlink(missing_ok=True)
            
            # Remove empty directories from deepest level up
            for dir_path in sorted(internal_dir.rglob('*'), reverse=True):
                if dir_path.is_dir():
                    dir_path.rmdir()
            
            # Remove the root internal directory
            internal_dir.rmdir()
            logging.info(f"Cleaned up internal working directory using alternative method: {internal_dir}")
            
        except Exception as e2:
            logging.error(f"Failed to clean up directory {internal_dir} using alternative method: {e2}")

def main():
    """Main execution function"""
    # Parse arguments
    # args are in ./cressent/modules/args.py
    args_ns = args.parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args_ns.outputDir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logger after creating output directory
    logger = setup_logger(output_dir)
    start_time = time.time()
    
    # Setup internal working folders
    internal_dir = output_dir / "internal"
    input_dir = internal_dir / "gff-split-inputs"
    internal_output_dir = internal_dir / "gff-split-outputs"
    
    # Create internal directories
    input_dir.mkdir(parents=True, exist_ok=True)
    internal_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Determine the input FASTA full path
        def validate_fasta(filename):
            with open(filename, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                if any(fasta):
                    print("FASTA checked.")
                    input_fasta = filename
                    return input_fasta
                else:
                    sys.exit("Error: Input file is not in the FASTA format.\n")
                    logging.info(f"Using: {input_fasta} is not in the FASTA format")
        # Convert input paths to absolute paths
        input_fasta = Path(args_ns.inputFasta).resolve()
        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
        # check fasta
        input_fasta = validate_fasta(input_fasta)
        input_gff = Path(args_ns.inputGFF).resolve()
        # Verify input files exist
        if not input_gff.exists():
            raise FileNotFoundError(f"Input GFF file not found: {input_gff}")
        output_gff = output_dir / Path(args_ns.outputGFF).name
                
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Input FASTA: {input_fasta}")
        logger.info(f"Input GFF: {input_gff}")
        logger.info(f"Output GFF: {output_gff}")
        logger.info(f"Internal working directory: {internal_dir}")
            
        # Process files
        input_convert(input_fasta, input_gff, input_dir)
        
        stats = {
            'total_files': 0,
            'files_with_iterons': 0,
            'total_iterons': 0
        }
        
        # Process each input file
        for entry in input_dir.glob('*.gff'):
            output_file = cruise.getOutputFile(entry, internal_output_dir)
            
            logger.info(f"Processing {entry.name}")
            result = cruise.findPosIterons(
                str(entry),  # Convert Path to string for compatibility
                output_file,
                args_ns.minLength,
                args_ns.maxLength,
                args_ns.range,
                args_ns.wiggle,
                args_ns.rank,
                args_ns.numberTopIterons,
                args_ns.maxScore,
                args_ns.goodLength,
                args_ns.doStemLoop,
                args_ns.doKnownIterons,
                args_ns.maxDist,
                args_ns.bestDist,
                args_ns.scoreRange
            )
            
            if result is not None:
                found_iterons, new_iteron_count = result
                stats['total_files'] += 1
                
                if found_iterons:
                    stats['files_with_iterons'] += 1
                    stats['total_iterons'] += new_iteron_count
        
        # Combine outputs
        output_convert(internal_output_dir, output_gff, args_ns.outputAnnotations)
        
        # Report results
        success_rate = (stats['files_with_iterons'] / stats['total_files'] * 100 
                       if stats['total_files'] > 0 else 0)
        
        logger.info(
            f"\nAnalysis complete:\n"
            f"- Found {stats['total_iterons']} iteron candidates\n"
            f"- In {stats['files_with_iterons']} out of {stats['total_files']} genomes\n"
            f"- Success rate: {success_rate:.1f}%\n"
            f"- Internal files preserved in: {internal_dir}"
        )
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise

    finally:
        # Clean up internal directories after processing is complete
        cleanup_internal_dirs(internal_dir)
        
    elapsed_time = time.time() - start_time
    logger.info(f"Total execution time: {elapsed_time:.3f} seconds")

if __name__ == "__main__":
    main()


