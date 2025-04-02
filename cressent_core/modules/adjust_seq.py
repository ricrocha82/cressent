import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def is_nucleotide_sequence(sequence: str) -> bool:
    """Check if a sequence contains only valid nucleotide characters (A, T, G, C, and optionally N)."""
    return bool(re.fullmatch(r'^[ATGCNatgcn]+$', sequence))

def adjust_sequence_start(sequence, motif="TAGTATTAC"):
    """Rotate the sequence so it starts with the given motif."""
    motif_index = sequence.find(motif)
    
    if motif_index == -1:
        print(f"Motif '{motif}' not found in sequence: {sequence[:30]}...")  # First 30 bases shown
        return sequence  # Return original if motif not found
    
    print(f"Motif '{motif}' found. Adjusting sequence start position.")
    return sequence[motif_index:] + sequence[:motif_index]

def process_fasta(input_fasta, output_dir, motif="TAGTATTAC"):
    """Adjust sequences in a FASTA file to start with the specified motif."""

    # Generate output file path
    prefix = os.path.splitext(os.path.basename(input_fasta))[0]
    output_fasta = os.path.join(output_dir, f"{prefix}_motif_adj.fa")

    with open(output_fasta, "w") as out_fh:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_str = str(record.seq)
            
            if not is_nucleotide_sequence(seq_str):
                print(f"Skipping non-nucleotide sequence: {record.id}")
                continue
            
            adjusted_seq = adjust_sequence_start(seq_str, motif)
            record.seq = Seq(adjusted_seq)  # Ensure sequence type compatibility
            SeqIO.write(record, out_fh, "fasta")

    return output_fasta  # Return the adjusted file path for downstream use

def main():
    parser = argparse.ArgumentParser(description="Adjust sequences in a FASTA file to start with a specified motif.")
    parser.add_argument("-i","--input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("-o","--output_fasta", default=".", help="Path to the output directory.")
    parser.add_argument("-m","--motif", default="TAGTATTAC", help="Motif to adjust sequences to start with (default: TAGTATTAC).")
    
    args = parser.parse_args()
    output_dir = args.output_fasta
    os.makedirs(output_dir , exist_ok=True)

    def validate_fasta(filename):
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                print("FASTA checked.")
                input_fasta = filename
                return input_fasta
            else:
                sys.exit("Error: Input file is not in the FASTA format.\n")

    # check fasta
    input_fasta = validate_fasta(args.input_fasta)
    process_fasta(input_fasta, output_dir , args.motif)

if __name__ == "__main__":
    main()



# python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/adjust_seq.py \
#             -i /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/rearrange/alpha_complete.fa \
#             -o /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/rearrange1 \
#               -m "CCGCAAATAACACTAAC"
