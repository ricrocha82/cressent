#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "plot_tree.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def main():
    # Define command-line arguments for the Python wrapper.
    parser = argparse.ArgumentParser(
        description="Wrapper to run the plot_tree_cli.R script from the ssDNA tool."
    )
    parser.add_argument("-t", "--tree", help="Input tree file (Newick format)")
    parser.add_argument("--dist_matrix", help="Use distance matrix method for tree construction")
    parser.add_argument("-o", "--output", required=True, default=".",help="Path to the output directory (Default: working directory)")
    parser.add_argument("--metadata_1", help="Optional CSV metadata file")
    parser.add_argument("--metadata_2", help="Optional TSV name table file")
    parser.add_argument("--alignment", help="Optional alignment file (FASTA) to include in the plot")
    parser.add_argument("--layout", default="rectangular",  help="Tree layout (e.g., rectangular, circular, unrooted) Default: rectangular")
    parser.add_argument("--branch_length", default="branch.length", help="Branch length parameter for ggtree (Default: branch.length)")
    parser.add_argument("--open_angle", default=0, type=float, help="Open angle for circular/unrooted layouts (default = 0)")
    parser.add_argument("--offset", default=0.14, type=float, help="Tip label offset (default = 0)")
    parser.add_argument("--tip_label", default="family", help="Column name to use as the tip label (default = family)")
    parser.add_argument("--color", default=True, type=lambda s: s.lower() not in ['false', '0', 'no'], 
                            help="Color tree by group (requires metadata) (default = True)")
    parser.add_argument("--fig_width", type=float, default=7, help="Figure width (default = 7)")
    parser.add_argument("--fig_height", type=float, default=7, help="Figure height (default = 7)")
    parser.add_argument("--plot_tips", default=True, type=lambda s: s.lower() not in ['false', '0', 'no'], 
                            help="Include tip labels in the plot (default = True)")
    parser.add_argument("--plot_name", default="tree_plot.pdf", help="Name of the output plot file")
    args = parser.parse_args()
    
    # Determine the absolute path of the modules directory.
    # Assuming the directory structure:
    # ssDNA_tool/
    #   ssDNA_annotator/
    #     modules/
    #       plot_tree.R
    modules_dir = Path(__file__).resolve().parent  # Directory where the script is located
    r_script_path = modules_dir / "plot_tree.R"  # Update with the actual path
    
    # Set up logging with log file in the specified output directory.
    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    setup_logging(args.output)
    logger = logging.getLogger("plot_tree.log")

    if not os.path.isfile(r_script_path):
        logger.error("R script 'plot_tree.R' not found in modules directory: %s", r_script_path)
        sys.exit(1)
    
    # Build a list of arguments to pass to the R script in the form --option=value
    r_args = []
    if args.tree is not None:
        r_args.append(f"--tree={args.tree}")
    if args.dist_matrix is not None:
        r_args.append(f"--dist_matrix={args.dist_matrix}")
    r_args.append(f"--outdir={args.output}")
    r_args.append(f"--layout={args.layout}")
    r_args.append(f"--branch_length={args.branch_length}")
    r_args.append(f"--open_angle={args.open_angle}")
    r_args.append(f"--offset={args.offset}")
    r_args.append(f"--tip_label={args.tip_label}")
    r_args.append(f"--fig_width={args.fig_width}")
    r_args.append(f"--fig_height={args.fig_height}")
    r_args.append(f"--plot_name={args.plot_name}")
    r_args.append(f"--color={str(args.color).upper()}")
    r_args.append(f"--plot_tips={str(args.plot_tips).upper()}")
    
    # Add optional arguments if provided
    if args.metadata_1:
        r_args.append(f"--metadata_1={args.metadata_1}")
    if args.metadata_2:
        r_args.append(f"--metadata_2={args.metadata_2}")
    if args.alignment:
        r_args.append(f"--alignment={args.alignment}")
    
    # Construct the full command: Rscript path/to/plot_tree_cli.R <arguments>
    command = ["Rscript", r_script_path] + r_args
    logger.info("Executing command:")
    logger.info(" ".join(map(str, command)))
    
    # Execute the R script with the specified arguments
    try:
        subprocess.run(command, check=True)
        logger.info("R script executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error("R script execution failed with error: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()



