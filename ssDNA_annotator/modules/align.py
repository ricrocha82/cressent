#!/usr/bin/env python3

import argparse
import subprocess
import os
import re

def run_command(command, error_message):
    """
    Helper function to run a shell command and handle errors.
    """
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command succeeded: {command}")
    except subprocess.CalledProcessError as e:
        print(f"{error_message}: {e}")
        exit(1)


def align_sequences(input_fasta, output_dir, prefix, mafft_ep=0.123, threads=1):
    """
    Align sequences using MAFFT and save output in the output directory.
    """

    print(f"\n\n Alignment using MAFFT \n\n")
    # aligned_file = os.path.join(output_dir, "aligned_sequences.fasta")
    aligned_file = os.path.join(output_dir, f"{prefix}_aligned_sequences.fasta")
    cmd = f"mafft --localpair --maxiterate 1000 --leavegappyregion --thread {threads} --ep {mafft_ep} {input_fasta} > {aligned_file}"
    run_command(cmd, "Error running MAFFT")
    return aligned_file


def trim_alignment(aligned_fasta, output_dir, prefix, gap_threshold=0.2):
    """
    Trim alignment using TrimAl.
    """
    print(f"\n\n Trimming alignment using trimal \n\n")

    # trimmed_file = os.path.join(output_dir, "aligned_trimmed_sequences.fasta")
    trimmed_file = os.path.join(output_dir, f"{prefix}_aligned_trimmed_sequences.fasta")
    cmd = f"trimal -keepheader -in {aligned_fasta} -out {trimmed_file} -gt {gap_threshold}"
    run_command(cmd, "Error running TrimAl")
    return trimmed_file

def main():
    parser = argparse.ArgumentParser(description="Pipeline for sequence alignment, trimming, motif splitting, and phylogenetic tree generation.")
    parser.add_argument("-t","--threads", default=1, help="Number of threads")
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing input FASTA file and for saving outputs.")
    # MAFFT args
    parser.add_argument("--mafft_ep", type=float, default=0.123, help="alignment length for MAFFT aligment (default: 0.123).")
    # TRIMAL args
    parser.add_argument("--gap_threshold", type=float, default=0.2, help="Gap threshold for TrimAl trimming (default: 0.2).")

    args = parser.parse_args()
    prefix = os.path.splitext(os.path.basename(args.input_fasta))[0]
    input_fasta = os.path.join(args.directory, args.input_fasta)
    output_dir = args.directory

    # Create nested directories
    try:
        os.makedirs(output_dir)
        print(f"directory '{output_dir}' created successfully.")
    except FileExistsError:
        print(f"directory '{output_dir}' already exist.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{nested_directory}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # Run the pipeline
    aligned_fasta = align_sequences(input_fasta, output_dir, prefix, mafft_ep=args.mafft_ep, threads=args.threads)
    trimmed_fasta = trim_alignment(aligned_fasta, output_dir, prefix, gap_threshold=args.gap_threshold)
    print(f"\n\n Filtered alignment is saved in '{output_dir} \n\n")


if __name__ == "__main__":
    main()


#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/reps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/

#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_test/test_data/Circoviridae_caps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_3/

#  /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py --threads 24 \
#                 --input_fasta /fs/project/PAS1117/ricardo/ssDNA_test/test_data/Circoviridae_reps.fa \
#                 -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_3/