#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate gene arrow plots from motif data using R')
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True, help='Input file path for the motif table CSV')
    
    # Optional arguments
    parser.add_argument('-o', '--output', default='gene_motif.pdf', help='Output file path for the generated plot (default: gene_motif.pdf)')
    parser.add_argument('--height', type=float, default=10, help='Height of the output plot in inches (default: 10)')
    parser.add_argument('--width', type=float, default=10, help='Width of the output plot in inches (default: 10)')
    parser.add_argument('-t', '--title', help='Title for the plot (optional)')
    
    # Parse arguments
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    
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
    cmd = ["Rscript", r_script_path, args.input, args.output, str(args.height), str(args.width)]
    
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
            for line in process.stderr:
                print(line.strip())
            sys.exit(process.returncode)
        
        print(f"Successfully created plot: {args.output}")
        
    except Exception as e:
        print(f"Error executing R script: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()


# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/gene_map.py \
#                     -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc/motif_table.csv \
#                     -o /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc/gene_map_motif_py.pdf \
#                     --height 10 --width 10 \
#                     -t "Motif Distribution in Genes"