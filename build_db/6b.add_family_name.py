# source ~/miniconda3/bin/activate
# conda activate biopython

import os
import glob
from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser(description='Update FASTA headers with family names')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing FASTA files')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for modified FASTA files')
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for fasta_file in glob.glob(os.path.join(args.input_dir, "*.fa")):
        family_name = os.path.basename(fasta_file).replace(".fa", "")
        output_file = os.path.join(args.output_dir, os.path.basename(fasta_file))
        
        with open(output_file, "w") as out_f:
            for record in SeqIO.parse(fasta_file, "fasta"):
                record.id = f"{record.description}|{family_name}"
                record.description = ""
                SeqIO.write(record, out_f, "fasta")

    print(f"Processing complete! Files saved in: {args.output_dir}")

if __name__ == "__main__":
    main()