#!/usr/bin/env python3

import argparse
import subprocess
import os
import re
from pathlib import Path


def run_command(command, error_message):
    """
    Helper function to run a shell command and handle errors.
    """
    try:
        subprocess.run(command, shell=False, check=True)
        print(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        print(f"{error_message}: {e}")
        exit(1)

def generate_tanglegram(tree1, tree2, label1, label2, output_dir, name_tanglegram="tanglegram.pdf", width=20, height=11,lab_cex=1.5):
    """
    Generate a tanglegram using an R script.
    
    Args:
        tree1 (str): Path to the first tree file.
        tree2 (str): Path to the second tree file.
        label1 (str): Label for the first tree.
        label2 (str): Label for the second tree.
        output_dir (str): Directory to save the output.
        name_tanglegram (str): Name of the output PDF file (default: "tanglegram.pdf").
        width (numeric): width of the pdf (default: 20 units).
        height (numeric): height of the pdf (default: 11 units).
        lab_cex (numeric): cex size of the labels.
    """
    # Path to the R script
    modules_dir = Path(__file__).resolve().parent  # Directory where the script is located
    r_script_path = modules_dir / "tanglegram.R"  # Update with the actual path
    #r_script_path = "/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/tanglegram.R"  # Update with the actual path

    # Build the command
    cmd = [
        "Rscript",
        r_script_path,
        f"tree1={tree1}",
        f"tree2={tree2}",
        f"label1={label1}",
        f"label2={label2}",
        f"output={output_dir}",
        f"name_tanglegram={name_tanglegram}",
        f"height={height}",
        f"width={width}"
        f"lab.cex={lab_cex}"
    ]

    # Run the R script
    run_command(cmd, "Error making tanglegram")

def main():
    parser = argparse.ArgumentParser(description = "Generate a tanglegram from two phylogenetic trees.")
    parser.add_argument("--tree1", required = True, help = "Path to the first tree file.") # arg_tree1
    parser.add_argument("--tree2", required = True, help = "Path to the second tree file.") # arg_tree2
    parser.add_argument("--label1", required = True, help = "Label for the first tree in the tanglegram.") # arg_label1
    parser.add_argument("--label2", required = True, help = "Label for the second tree in the tanglegram.") # arg_label2
    parser.add_argument("--output", required = True, help = "Directory where the tanglegram will be saved.") # arg_output
    parser.add_argument("--name_tanglegram", default = "tanglegram.pdf", help = "Name of the tanglegram PDF file.") 
    parser.add_argument("--width", default = 20, type=float, help = "Width of the tanglegram PDF file.") 
    parser.add_argument("--height", default = 11, type=float, help = "Height of the tanglegram PDF file.") 
    parser.add_argument("--lab_cex", default = 1.5, type=float, help = "cex size of the labels.") 

    args = parser.parse_args()
    
    # run R script
    generate_tanglegram(args.tree1, args.tree2, args.label1, args.label2, args.output, name_tanglegram=args.name_tanglegram, width=args.width, height=args.height, lab_cex=args.lab_cex)


if __name__ == "__main__":
    main()


# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/tanglegram.py \
#     --tree1 /fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree_helic/sanitized_sequences.fasta.treefile \
#     --tree2 /fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree_endonuc/sanitized_sequences.fasta.treefile \
#     --label1 Helicase \
#     --label2 Endonuclease \
#     --output /fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/ \
#     --name_tanglegram "my_tanglegram.pdf" \
#     --height 11 \
#     --width 20 \
#     --lab_cex 2