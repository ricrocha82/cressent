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

def find_pattern(input_fasta, seq_pattern, output_dir,table_name="pattern_positions.txt"):
    """
    Find positions based on the pattern using seqkit.
    """
    seq_table = os.path.join(output_dir, table_name)

    # Find pattern positions using seqkit
    cmd = f"seqkit locate --id-ncbi -i -r -p '{seq_pattern}' {input_fasta} > {seq_table}"
    run_command(cmd, "Error locating patterns with SeqKit")

def main():
    parser = argparse.ArgumentParser(description="Pipeline to split a sequence based on determined patter (e.g., motif)")
    parser.add_argument("-i","--input_fasta", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("-d", "--directory", required=True, help="Directory for saving outputs.")
    parser.add_argument("-p", "--pattern", required=True, help="Sequence pattern (regex) for splitting sequences.")
    parser.add_argument("-n","--table_name", default="pattern_positions.txt", help="Sequence pattern (regex) for splitting sequences.")

    args = parser.parse_args()
    input_fasta = os.path.join(args.directory, args.input_fasta)
    output_dir = args.directory

    # Create nested directories
    try:
        os.makedirs(output_dir)
        print(f"directory '{output_dir}' created successfully.")
    except FileExistsError:
        print(f"directory '{output_dir}' already exist.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{output_dir}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # run pipeline
    find_pattern(input_fasta, args.pattern, output_dir, table_name=args.table_name)

    print(f"Pattern found. Outputs saved in {output_dir}")


if __name__ == "__main__":
    main()



# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/find_pattern.py \
#             -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/aligned_trimmed_sequences.fasta \
#             --pattern [GA].{4}GK[TS] \
#             -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/ \
#             -n position_text_2.txt

# /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/find_pattern.py \
#             -i /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/families_nc/Anelloviridae/Anelloviridae_rep.fasta \
#             --pattern AAGTATT*AC \
#             -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/ \
#             -n Anelloviridae_rep_position.txt

