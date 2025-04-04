#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate gene arrow plots from motif data using R')

# Required arguments
parser.add_argument('-i', '--input', required=True, help='Input file path for the motif table CSV')

# Optional arguments
parser.add_argument('-o', '--output', default='./', help='Output directory where the plot will be saved')
parser.add_argument('--filename', default='gene_motif.pdf', help='Name of the output file (default: gene_motif.pdf)')
parser.add_argument('--height', type=float, default=10, help='Height of the output plot in inches (default: 10)')
parser.add_argument('--width', type=float, default=10, help='Width of the output plot in inches (default: 10)')
parser.add_argument('-t', '--title', help='Title for the plot (optional)')

def main():
    # Parse arguments
    args = parser.parse_args()
    
    # Ensure filename has a proper extension
    if not os.path.splitext(args.filename)[1]:
        args.filename = args.filename + '.pdf'
        print(f"No file extension detected in filename, using: {args.filename}")
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(args.output):
        print(f"Creating output directory: {args.output}")
        os.makedirs(args.output, exist_ok=True)
    
    # Combine output directory and filename
    output_file_path = os.path.join(args.output, args.filename)
    print(f"Output will be saved to: {output_file_path}")
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file does not exist: {args.input}")
        sys.exit(1)
    
    # Get the path to the R script - assuming it's in the same directory as this script
    r_script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_map.R")
    
    if not os.path.exists(r_script_path):
        print(f"Error: R script not found at {r_script_path}")
        sys.exit(1)
    
    # Prepare command
    cmd = ["Rscript", r_script_path, args.input, output_file_path, str(args.height), str(args.width)]
    
    # Add title if provided
    if args.title:
        cmd.append(args.title)
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        # Run the R script
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Stream output in real-time
        for line in process.stdout:
            print(line.strip())
        
        # Wait for the process to complete
        process.wait()
        
        # Check if process was successful
        if process.returncode != 0:
            print("Error: R script execution failed")
            stderr_output = process.stderr.read()
            print(stderr_output)
            sys.exit(process.returncode)
        
        print(f"Successfully created plot: {output_file_path}")
        
    except Exception as e:
        print(f"Error executing R script: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()