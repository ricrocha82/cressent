#!/usr/bin/env python3

"""
Recombination detection module for CRESSENT.

This module integrates OpenRDP (https://github.com/aglucaci/OpenRDP) for detecting recombination events in ssDNA virus sequences.

Copyright (C) 2025 CRESSENT developers

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
import shutil
import tempfile
import logging
import zipfile
import tarfile
import platform
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
import sys
import importlib.util

# Get the absolute path to the module directory
# Get the absolute path to the module directory
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

# Check if we're in the expected directory structure
if "modules" in MODULE_DIR and os.path.exists(os.path.join(MODULE_DIR, "openrdp")):
    OPENRDP_DIR = os.path.join(MODULE_DIR, "openrdp")
elif os.path.exists(os.path.join(MODULE_DIR, "modules", "openrdp")):
    OPENRDP_DIR = os.path.join(MODULE_DIR, "modules", "openrdp")
else:
    # Try to find openrdp directory
    OPENRDP_DIR = None
    for root, dirs, files in os.walk(MODULE_DIR):
        if "openrdp" in dirs:
            OPENRDP_DIR = os.path.join(root, "openrdp")
            break
    
    # If still not found, use default
    if OPENRDP_DIR is None:
        OPENRDP_DIR = os.path.join(MODULE_DIR, "openrdp")
        os.makedirs(OPENRDP_DIR, exist_ok=True)

def setup_logger(output_dir, log_level=logging.INFO, quiet=False):
    """
    Setup logger for recombination module.
    
    Parameters:
    -----------
    output_dir : str
        Directory where log file will be saved
    log_level : int
        Logging level (default: INFO)
    quiet : bool
        Whether to suppress console output
    
    Returns:
    --------
    logger : logging.Logger
        Configured logger
    """
    # Create logger
    logger = logging.getLogger('ssDNA_recombination')
    logger.setLevel(log_level)
    logger.handlers = []  # Clear any existing handlers
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Create file handler
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(output_dir, f'recombination_{timestamp}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # Create console handler if not quiet
    if not quiet:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    return logger

def check_binary_compatibility(binary_path, logger):
    """
    Check if a binary executable is compatible with the current system.
    
    Parameters:
    -----------
    binary_path : str
        Path to the binary executable
    logger : logging.Logger
        Logger for tracking the process
    
    Returns:
    --------
    bool
        True if compatible, False otherwise
    """
    if not os.path.exists(binary_path):
        logger.warning(f"Binary not found: {binary_path}")
        return False
        
    try:
        # Try to run the binary with --help or --version
        result = subprocess.run(
            [binary_path, "--help"], 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        # Check if the binary ran without GLIBC errors
        stderr = result.stderr.decode().lower()
        if 'glibc' in stderr:
            logger.warning(f"Binary compatibility issue detected for {binary_path}: {stderr}")
            return False
        return True
    except Exception as e:
        logger.warning(f"Binary compatibility check failed for {binary_path}: {e}")
        return False

def compile_3seq(source_zip, target_dir, logger):
    """
    Compile 3Seq from source code.
    
    Parameters:
    -----------
    source_zip : str
        Path to 3Seq source code zip file
    target_dir : str
        Directory where compiled binary should be placed
    logger : logging.Logger
        Logger for tracking the process
    
    Returns:
    --------
    bool
        True if compilation successful, False otherwise
    """
    temp_dir = tempfile.mkdtemp()
    logger.info(f"Extracting 3Seq source code to {temp_dir}")
    
    try:
        # Extract source code
        with zipfile.ZipFile(source_zip, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
        
        # Find the source directory
        dirs = [d for d in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, d))]
        if not dirs:
            source_dir = temp_dir
        else:
            source_dir = os.path.join(temp_dir, dirs[0])
        
        logger.info(f"Compiling 3Seq in {source_dir}")
        
        # Compile
        result = subprocess.run(
            ["make"],
            cwd=source_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Failed to compile 3Seq: {result.stderr}")
            return False
        
        # Find the compiled binary
        binary_paths = []
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if file == '3seq' or file.startswith('3seq.'):
                    binary_paths.append(os.path.join(root, file))
        
        if not binary_paths:
            logger.error("Compiled 3Seq binary not found")
            return False
        
        # Copy binary to target directory
        binary_path = binary_paths[0]
        target_binary = os.path.join(target_dir, "3seq.Unix")
        shutil.copy2(binary_path, target_binary)
        os.chmod(target_binary, 0o755)
        logger.info(f"Successfully compiled and installed 3Seq to {target_binary}")

        # Check if P-value table exists
        p_table_path = os.path.join(target_dir, "myPvalueTable")
        if not os.path.exists(p_table_path):
            logger.info("P-value table not found. Generating a new 500 × 500 × 500 table...")
            
            # Generate the P-value table
            try:
                # Change to target directory
                original_dir = os.getcwd()
                os.chdir(target_dir)
                
                # Run 3seq to generate the table
                gen_result = subprocess.run(
                    ["./3seq.Unix", "-g", "myPvalueTable", "500"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                if gen_result.returncode != 0:
                    logger.error(f"Failed to generate P-value table: {gen_result.stderr}")
                    logger.warning("3Seq may not work correctly without a P-value table")
                else:
                    logger.info("Successfully generated P-value table")
                
                # Return to original directory
                os.chdir(original_dir)

            except Exception as e:
                    logger.error(f"Error generating P-value table: {e}")
                    logger.warning("3Seq may not work correctly without a P-value table")
        else:
            logger.info(f"P-value table already exists at {p_table_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error compiling 3Seq: {e}", exc_info=True)
        return False
    
    finally:
        # Clean up
        shutil.rmtree(temp_dir)

def compile_geneconv(source_tar, target_dir, logger):
    """
    Compile GENECONV from source code using the specified compilation command.
    
    Parameters:
    -----------
    source_tar : str
        Path to GENECONV source code tar file
    target_dir : str
        Directory where compiled binary should be placed
    logger : logging.Logger
        Logger for tracking the process
    
    Returns:
    --------
    bool
        True if compilation successful, False otherwise
    """
    temp_dir = tempfile.mkdtemp()
    logger.info(f"Extracting GENECONV source code to {temp_dir}")
    
    try:
        # Extract source code
        if source_tar.endswith('.tar.gz'):
            with tarfile.open(source_tar, 'r:gz') as tar_ref:
                tar_ref.extractall(temp_dir)
        else:
            with tarfile.open(source_tar, 'r') as tar_ref:
                tar_ref.extractall(temp_dir)
        
        # Find the source directory
        source_dir = None
        for root, dirs, files in os.walk(temp_dir):
            # Look for key source files
            if 'geneconv.c' in files or 'GENECONV.C' in files:
                source_dir = root
                break
        
        if not source_dir:
            logger.error("Could not find GENECONV source directory in the tar file")
            return False
        
        logger.info(f"Compiling GENECONV in {source_dir}")
        
        # Check for required source files (case-insensitive)
        required_files = ['geneconv.c', 'version.c', 'vcalc.c', 'vtcalc.c', 
                         'vsetopts.c', 'vread.c', 'vdump.c', 'vutil.c']
        
        actual_files = {}  # Map lowercase names to actual filenames
        for file in os.listdir(source_dir):
            actual_files[file.lower()] = file
        
        source_files = []
        for req_file in required_files:
            if req_file.lower() in actual_files:
                source_files.append(actual_files[req_file.lower()])
            else:
                logger.error(f"Required source file {req_file} not found")
                return False
        
        # Run the specified compilation command
        logger.info("Running gcc compilation command")
        compile_cmd = ["gcc", "-DUNIX", "-o", "geneconv", "-O3"]
        compile_cmd.extend(source_files)
        compile_cmd.append("-lm")
        
        result = subprocess.run(
            compile_cmd,
            cwd=source_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode != 0:
            logger.warning(f"Compilation with -lm failed: {result.stderr}")
            
            # Try without math library
            compile_cmd.pop()  # Remove -lm
            logger.info("Trying compilation without -lm")
            result = subprocess.run(
                compile_cmd,
                cwd=source_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            if result.returncode != 0:
                logger.error(f"Compilation without -lm also failed: {result.stderr}")
                return False
        
        # Check if geneconv executable was created
        geneconv_path = os.path.join(source_dir, 'geneconv')
        if not os.path.exists(geneconv_path):
            logger.error("geneconv binary not found after compilation")
            return False
        
        # Copy binary to target directory
        target_binary = os.path.join(target_dir, "geneconv.Unix")
        shutil.copy2(geneconv_path, target_binary)
        os.chmod(target_binary, 0o755)
        logger.info(f"Successfully compiled and installed GENECONV to {target_binary}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error compiling GENECONV: {e}", exc_info=True)
        return False
    
    finally:
        # Clean up
        shutil.rmtree(temp_dir)

def ensure_binaries_exist(logger):
    """
    Check if binary executables exist, compile from source if not.
    
    Parameters:
    -----------
    logger : logging.Logger
        Logger for tracking the process
    """
    logger.info("Checking for required binary executables...")
    
    # Create binary directories if they don't exist
    os.makedirs(os.path.join(OPENRDP_DIR, "bin", "3Seq"), exist_ok=True)
    os.makedirs(os.path.join(OPENRDP_DIR, "bin", "GENECONV"), exist_ok=True)
    
    # Check platform
    system = platform.system()
    
    # Define expected binary files based on platform
    if system == "Windows":
        threeseq_bin = os.path.join(OPENRDP_DIR, "bin", "3Seq", "windows_3seq.exe")
        geneconv_bin = os.path.join(OPENRDP_DIR, "bin", "GENECONV", "windows_geneconv.exe")
    elif system == "Darwin":  # macOS
        threeseq_bin = os.path.join(OPENRDP_DIR, "bin", "3Seq", "3seq.macOS")
        geneconv_bin = os.path.join(OPENRDP_DIR, "bin", "GENECONV", "geneconv.macOS")
    else:  # Linux/Unix
        threeseq_bin = os.path.join(OPENRDP_DIR, "bin", "3Seq", "3seq.Unix")
        geneconv_bin = os.path.join(OPENRDP_DIR, "bin", "GENECONV", "geneconv.Unix")
    
    # Check if binaries already exist
    threeseq_exists = os.path.exists(threeseq_bin)
    geneconv_exists = os.path.exists(geneconv_bin)
    
    # Set paths to source code archives
    threeseq_source = os.path.join(OPENRDP_DIR, "bin", "source_code", "3seq_build_170612.zip")
    geneconv_source = os.path.join(OPENRDP_DIR, "bin", "source_code", "geneconv_unix.source.tar")
    
    # Check for compatibility issues
    if threeseq_exists:
        logger.info(f"3Seq binary found at {threeseq_bin}")
        if system != "Windows":
            try:
                os.chmod(threeseq_bin, 0o755)
            except Exception as e:
                logger.warning(f"Could not set executable permissions for 3Seq: {e}")
        
        # Check if binary is compatible with the system
        if not check_binary_compatibility(threeseq_bin, logger):
            logger.warning("3Seq binary is not compatible with your system. Will try to compile from source.")
            threeseq_exists = False
    
    if geneconv_exists:
        logger.info(f"GENECONV binary found at {geneconv_bin}")
        if system != "Windows":
            try:
                os.chmod(geneconv_bin, 0o755)
            except Exception as e:
                logger.warning(f"Could not set executable permissions for GENECONV: {e}")
        
        # Check if binary is compatible with the system
        if not check_binary_compatibility(geneconv_bin, logger):
            logger.warning("GENECONV binary is not compatible with your system. Will try to compile from source.")
            geneconv_exists = False
    
    # Compile binaries if needed
    if not threeseq_exists and os.path.exists(threeseq_source):
        logger.info(f"Compiling 3Seq from source: {threeseq_source}")
        threeseq_compiled = compile_3seq(
            threeseq_source,
            os.path.join(OPENRDP_DIR, "bin", "3Seq"),
            logger
        )
        if threeseq_compiled:
            logger.info("Successfully compiled 3Seq from source")
        else:
            logger.warning("Failed to compile 3Seq from source. 3Seq analysis will be disabled.")
    
    if not geneconv_exists and os.path.exists(geneconv_source):
        logger.info(f"Compiling GENECONV from source: {geneconv_source}")
        geneconv_compiled = compile_geneconv(
            geneconv_source,
            os.path.join(OPENRDP_DIR, "bin", "GENECONV"),
            logger
        )
        if geneconv_compiled:
            logger.info("Successfully compiled GENECONV from source")
        else:
            logger.warning("Failed to compile GENECONV from source. GENECONV analysis will be disabled.")
    
    # Final check
    threeseq_exists = os.path.exists(threeseq_bin) and check_binary_compatibility(threeseq_bin, logger)
    geneconv_exists = os.path.exists(geneconv_bin) and check_binary_compatibility(geneconv_bin, logger)
    
    return threeseq_exists, geneconv_exists

def get_platform_executable(tool_name):
    """
    Get the correct executable based on platform.
    
    Parameters:
    -----------
    tool_name : str
        Name of the tool ('3Seq' or 'GENECONV')
    
    Returns:
    --------
    str
        Path to the executable for the current platform
    """
    system = platform.system()
    bin_dir = os.path.join(OPENRDP_DIR, "bin")
    
    if tool_name == "3Seq":
        tool_dir = os.path.join(bin_dir, "3Seq")
        if system == "Windows":
            return os.path.join(tool_dir, "windows_3seq.exe")
        elif system == "Darwin":  # macOS
            return os.path.join(tool_dir, "3seq.macOS")
        else:  # Linux/Unix
            return os.path.join(tool_dir, "3seq.Unix")
    
    elif tool_name == "GENECONV":
        tool_dir = os.path.join(bin_dir, "GENECONV")
        if system == "Windows":
            return os.path.join(tool_dir, "windows_geneconv.exe")
        elif system == "Darwin":  # macOS
            return os.path.join(tool_dir, "geneconv.macOS")
        else:  # Linux/Unix
            return os.path.join(tool_dir, "geneconv.Unix")
    
    return None

def parse_args(args=None):
    """Parse command line arguments for recombination detection."""
    parser = argparse.ArgumentParser(description="Detect recombination events in ssDNA virus sequences using OpenRDP")
    
    parser.add_argument('-i', '--input', required=True, help='Input alignment file in FASTA format')
    parser.add_argument('-o', '--output', default='.', help='Path to the output directory (Default: working directory)')
    parser.add_argument('-f', '--output_file', required=True, help='Output file for results (CSV format)')
    parser.add_argument('-c', '--config', help='Configuration file in INI format for OpenRDP parameters')
    
    # Method options
    parser.add_argument('-rdp', action='store_true', help='Run RDP method')
    parser.add_argument('-threeseq', action='store_true', help='Run 3Seq method')
    parser.add_argument('-geneconv', action='store_true', help='Run GENECONV method')
    parser.add_argument('-maxchi', action='store_true', help='Run MaxChi method')
    parser.add_argument('-chimaera', action='store_true', help='Run Chimaera method')
    parser.add_argument('-bootscan', action='store_true', help='Run Bootscan method')
    parser.add_argument('-siscan', action='store_true', help='Run Siscan method')
    parser.add_argument('-all', action='store_true', help='Run all methods')
    
    parser.add_argument('-quiet', action='store_true', help='Suppress console output')
    parser.add_argument('-verbose', action='store_true', help='Enable verbose logging')
    
    return parser.parse_args(args)

def detect_recombination(input_file, output_file, output_dir='.', config_file=None, 
                         run_geneconv=False, run_three_seq=False, run_rdp=False,
                         run_siscan=False, run_maxchi=False, run_chimaera=False, 
                         run_bootscan=False, quiet=False, verbose=False):
    """
    Detect recombination in DNA sequences using OpenRDP.
    
    Parameters:
    -----------
    input_file : str
        Path to input alignment file in FASTA format
    output_file : str
        Path to output file for results (CSV format)
    output_dir : str
        Directory where all output files will be saved (default: current directory)
    config_file : str, optional
        Path to configuration file in INI format
    run_geneconv, run_three_seq, run_rdp, run_siscan, run_maxchi, run_chimaera, run_bootscan : bool
        Flag to run specific methods
    quiet : bool, optional
        Suppress console output
    verbose : bool, optional
        Enable verbose logging
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup logging
    log_level = logging.DEBUG if verbose else logging.INFO
    logger = setup_logger(output_dir, log_level, quiet)
    
    # Check for binary executables and compile if needed
    threeseq_available, geneconv_available = ensure_binaries_exist(logger)
    
    # Disable methods with unavailable binaries
    if run_three_seq and not threeseq_available:
        logger.warning("3Seq binary is not available. Disabling 3Seq analysis.")
        run_three_seq = False
    
    if run_geneconv and not geneconv_available:
        logger.warning("GENECONV binary is not available. Disabling GENECONV analysis.")
        run_geneconv = False
    
    # Import OpenRDP's main function
    logger.info("Importing OpenRDP modules...")
    
    try:
        # Ensure all parent directories are in the Python path
        module_parent_dir = os.path.dirname(MODULE_DIR)
        if module_parent_dir not in sys.path:
            sys.path.insert(0, module_parent_dir)
        
        if MODULE_DIR not in sys.path:
            sys.path.insert(0, MODULE_DIR)
            
        if OPENRDP_DIR not in sys.path:
            sys.path.insert(0, OPENRDP_DIR)
        
        # Ensure the scripts directory has an __init__.py file
        scripts_dir = os.path.join(OPENRDP_DIR, "scripts")
        os.makedirs(scripts_dir, exist_ok=True)
        init_file = os.path.join(scripts_dir, "__init__.py")
        if not os.path.exists(init_file):
            with open(init_file, "w") as f:
                pass  # Create an empty file
        
        # Create __init__.py for OPENRDP_DIR if it doesn't exist
        init_file = os.path.join(OPENRDP_DIR, "__init__.py")
        if not os.path.exists(init_file):
            with open(init_file, "w") as f:
                pass
        
        # First, try to import from an installed package
        try:
            # Try importing main directly
            from openrdp.main import openrdp
            logger.info("Successfully imported OpenRDP main function from installed package")
        except ImportError:
            # Try relative import within the package
            try:
                sys.path.append(os.path.join(OPENRDP_DIR))
                from main import openrdp
                logger.info("Successfully imported OpenRDP main function using direct import")
            except ImportError:
                # If relative import fails, try loading the module directly
                logger.info("Direct import failed, trying to load module directly")
                
                main_path = os.path.join(OPENRDP_DIR, "main.py")
                if not os.path.exists(main_path):
                    logger.error(f"Main module file not found at {main_path}")
                    # Look for main.py in other locations
                    for root, dirs, files in os.walk(OPENRDP_DIR):
                        if "main.py" in files:
                            main_path = os.path.join(root, "main.py")
                            logger.info(f"Found main.py at: {main_path}")
                            break
                
                if os.path.exists(main_path):
                    import importlib.util
                    spec = importlib.util.spec_from_file_location("openrdp.main", main_path)
                    main_module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(main_module)
                    openrdp = main_module.openrdp
                    logger.info("Successfully imported OpenRDP main function using spec loader")
                else:
                    raise ImportError(f"Cannot find main.py in OpenRDP directory: {OPENRDP_DIR}")
        
        # Log analysis parameters
        logger.info(f"Starting recombination analysis on {input_file}")
        logger.info(f"Output file: {output_file}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Config file: {config_file if config_file else 'Default configuration'}")
        
        # Log selected methods
        methods = []
        if run_geneconv:
            methods.append("GENECONV")
        if run_three_seq:
            methods.append("3Seq")
        if run_rdp:
            methods.append("RDP")
        if run_siscan:
            methods.append("Siscan")
        if run_maxchi:
            methods.append("MaxChi")
        if run_chimaera:
            methods.append("Chimaera")
        if run_bootscan:
            methods.append("Bootscan")
        
        if not methods:
            logger.info("No specific methods selected, running all available methods")
        else:
            logger.info(f"Selected methods: {', '.join(methods)}")
        
        # Make output_file path absolute if it's not already
        if not os.path.isabs(output_file):
            output_file = os.path.join(output_dir, output_file)
        
        # Change current working directory to output_dir for OpenRDP to save any additional files there
        original_dir = os.getcwd()
        os.chdir(output_dir)
        
        try:
            # Run OpenRDP
            logger.info("Running OpenRDP analysis...")
            openrdp(input_file, output_file, config_file, 
                    run_geneconv, run_three_seq, run_rdp,
                    run_siscan, run_maxchi, run_chimaera, 
                    run_bootscan, quiet)
            
            # Read and log results summary
            try:
                results_df = pd.read_csv(output_file)
                num_events = len(results_df)
                logger.info(f"Analysis complete. Found {num_events} potential recombination events.")
                
                # Log method summary
                if num_events > 0:
                    method_counts = results_df['Method'].value_counts()
                    for method, count in method_counts.items():
                        logger.info(f"  - {method}: {count} events")
            except Exception as e:
                logger.warning(f"Could not read results file for summary: {e}")
            
            logger.info(f"Results saved to {output_file}")
            return True
            
        finally:
            # Restore original working directory
            os.chdir(original_dir)
            
    except Exception as e:
        logger.error(f"Error during recombination detection: {e}", exc_info=verbose)
        return False

def main(args=None):
    """Main function to run recombination detection."""
    parsed_args = parse_args(args)
    
    # Check if input file exists
    if not os.path.exists(parsed_args.input):
        print(f"Error: Input file {parsed_args.input} not found.")
        return 1
    
    # Create output directory if it doesn't exist
    os.makedirs(parsed_args.output, exist_ok=True)
    
    # Determine methods to run
    run_geneconv = parsed_args.geneconv
    run_three_seq = parsed_args.threeseq
    run_rdp = parsed_args.rdp
    run_siscan = parsed_args.siscan
    run_maxchi = parsed_args.maxchi
    run_chimaera = parsed_args.chimaera
    run_bootscan = parsed_args.bootscan
    
    # If all is selected, run all methods
    if parsed_args.all:
        run_geneconv = True
        run_three_seq = True
        run_rdp = True
        run_siscan = True
        run_maxchi = True
        run_chimaera = True
        run_bootscan = True
    
    # If no methods selected, run all
    if not (run_geneconv or run_three_seq or run_rdp or run_siscan or 
            run_maxchi or run_chimaera or run_bootscan):
        run_geneconv = True
        run_three_seq = True
        run_rdp = True
        run_siscan = True
        run_maxchi = True
        run_chimaera = True
        run_bootscan = True

    output_dir = parsed_args.output
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
    # check fasta
    input_fasta = validate_fasta(parsed_args.input)
    
    # Run recombination detection
    success = detect_recombination(
        input_fasta,
        parsed_args.output_file,
        parsed_args.output,
        parsed_args.config,
        run_geneconv,
        run_three_seq,
        run_rdp,
        run_siscan,
        run_maxchi,
        run_chimaera,
        run_bootscan,
        parsed_args.quiet,
        parsed_args.verbose
    )
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())


