#!/usr/bin/env python3

# this is a modified version of StemLoop-Finder from 
# Pratt, A.A., Torrance, E.L., Kasun, G.W., Stedman, K.M. and de la Higuera, I., 2021. StemLoop-Finder: a tool for the detection of DNA hairpins with conserved motifs. Microbiology Resource Announcements, 10(26), pp.10-1128.
# https://doi.org/10.1128/mra.00424-21
# CIte if using this module
# All rights reserved.

import argparse
import csv
import logging
import os
import sys
from pathlib import Path
from Bio import SeqIO
import sys

# Add the parent directory to sys.path to help find pos_stem_loop
current_dir = Path(__file__).resolve().parent
if str(current_dir) not in sys.path:
    sys.path.append(str(current_dir))

try:
    import RNA
except ImportError:
    print("Error: ViennaRNA package (RNA) not installed. Please install with: pip install ViennaRNA")
    sys.exit(1)

try:
    import gffutils
except ImportError:
    print("Error: gffutils not installed. Please install with: pip install gffutils")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print("Error: Biopython not installed. Please install with: pip install biopython")
    sys.exit(1)

try:
    from pos_stem_loop import PosStemLoop
except ImportError:
    print("Error: Cannot find pos_stem_loop.py. Make sure it's in the same directory as this script.")
    sys.exit(1)

def setup_logger(output_dir):
    """Set up a logger that writes to a file in the output directory."""
    logger = logging.getLogger('sl_finder')
    logger.setLevel(logging.DEBUG)
    log_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Define log file path in the output directory
    log_file = os.path.join(output_dir, 'sl_finder.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(log_formatter)
    logger.addHandler(file_handler)
    return logger

# Helper functions (adapted from your original code)
# Modified IUPAC dictionary to handle both cases
iupac_to_nuc = {
    'a': 'a', 't': 'tu', 'g': 'g', 'c': 'c',
    'r': 'ag', 'y': 'ctu', 's': 'gc', 'w': 'atu',
    'k': 'gtu', 'm': 'ac', 'b': 'cgtu', 'd': 'agtu',
    'h': 'actu', 'v': 'acg', 'n': 'acgtu', 'u': 'ut',
    'A': 'a', 'T': 'tu', 'G': 'g', 'C': 'c',
    'R': 'ag', 'Y': 'ctu', 'S': 'gc', 'W': 'atu',
    'K': 'gtu', 'M': 'ac', 'B': 'cgtu', 'D': 'agtu',
    'H': 'actu', 'V': 'acg', 'N': 'acgtu', 'U': 'ut'
}

def reverse_complement(bases, isDNA=True):
    """Handle both upper and lower case bases."""
    if isDNA:
        complement = {
            'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm',
            'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
            'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b',
            'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
            'n': 'n', 'N': 'N', 'w': 'w', 'W': 'W'
        }
    else:
        complement = {
            'a': 'u', 'c': 'g', 'g': 'c', 'u': 'a',
            'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',
            'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm',
            'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
            'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b',
            'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
            'n': 'n', 'N': 'N'
        }
    return ''.join([complement.get(base, base) for base in bases[::-1]])

def ambiguitiesToBases(ambiguousDNA):
    return [iupac_to_nuc[base] for base in ambiguousDNA]

def equivalentTo(ambiguousDNA, givenSeq):
    """Compare sequences case-insensitively."""
    for i in range(len(givenSeq)):
        if givenSeq[i].lower() not in iupac_to_nuc[ambiguousDNA[i].lower()]:
            return False
    return True

def removeDuplicates(pos_loops):
    new_pos_loops = []
    if pos_loops:
        new_pos_loops.append(pos_loops[0])
    for loop in pos_loops[1:]:
        if loop != new_pos_loops[-1]:
            new_pos_loops.append(loop)
    return new_pos_loops

def bestNona(family_type):
    mapping = {
        "geminiviridae": "trakattrc",
        "circoviridae": "tagtattac",
        "cycloviridae": "taatattac",
        "general": "tagtattac",
        "genomoviridae": "tawwdhwan",
        "smacoviridae": "nakwrttac",
    }
    return mapping.get(family_type.lower(), "nantantan")

def motifSearch(sequence, name, motif, basesAroundMotif):
    """Search for motifs in sequence, handling both cases."""
    pos_motifs = []
    motif_len = len(motif)
    # Convert motif to lowercase for consistent comparison
    motif = motif.lower()
    
    for i in range(basesAroundMotif, len(sequence) - (motif_len - 1)):
        subseq = sequence[i:i + motif_len]
        expanded_start = i - basesAroundMotif
        expanded_end = i + basesAroundMotif + motif_len
        expanded_subseq = sequence[expanded_start:expanded_end]
        
        # Store original case for output
        original_subseq = subseq
        # Compare in lowercase
        if equivalentTo(motif, subseq.lower()):
            pos_motifs.append(PosStemLoop(name, expanded_subseq, expanded_start, original_subseq, i + 1, False))
        elif equivalentTo(reverse_complement(motif, True), subseq.lower()):
            pos_motifs.append(PosStemLoop(name, expanded_subseq, expanded_start, reverse_complement(original_subseq, True), i + 1, True))
    return pos_motifs

def nonaSearch(sequence, name, family, basesAroundNona):
    if family:
        return motifSearch(sequence, name, bestNona(family), basesAroundNona)
    else:
        return motifSearch(sequence, name, "nantantan", basesAroundNona)

def run_sl_finder(fasta_file, gff_file, out_gff, output_dir, out_csv=None, **kwargs):
    """
    Run the SL-Finder module.
    """
    # Validate input files exist
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_file}")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"Input GFF file not found: {gff_file}")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up logger
    logger = setup_logger(output_dir)
    logger.info(f"Starting SL-Finder module with input files: {fasta_file}, {gff_file}")
    
    # Parse parameters
    motif = kwargs.get('motif', "nantantan")
    family = kwargs.get('family', None)
    ideal_stem_len = kwargs.get('idealstemlen', 11)
    ideal_loop_len = kwargs.get('ideallooplen', 11)
    frame = kwargs.get('frame', 15)
    
    try:
        # Read FASTA using Biopython
        fasta_records = list(SeqIO.parse(fasta_file, "fasta"))
        if not fasta_records:
            raise ValueError("No sequences found in FASTA file")
        logger.info(f"Parsed {len(fasta_records)} FASTA records.")
        
        # Create temporary GFF database
        db = gffutils.create_db(gff_file, ":memory:", merge_strategy='create_unique', 
                              keep_order=True, force=True)
        
        features_out = []
        csv_rows = []
        
        for record in fasta_records:
            seq_str = str(record.seq).lower() + str(record.seq)[:50]
            logger.info(f"Processing record: {record.id}")
            
            # Search for stem loops
            if family:
                pos_loops = nonaSearch(seq_str, record.id, family, frame)
            else:
                pos_loops = motifSearch(seq_str, record.id, motif, frame)
            
            # Process each potential stem loop
            for pos_loop in pos_loops:
                try:
                    folded = pos_loop.fold()
                    stem_loop_seq = pos_loop.find_stem_loop(ideal_stem_len, ideal_loop_len)
                    
                    if family:
                        pos_loop.diff_score_modify(bestNona(family))
                    
                    if (pos_loop.stem_start != -1 and 
                        pos_loop.score < 15 and 
                        pos_loop.stem_end > (pos_loop.motif_start + 12) and 
                        pos_loop.motif_start > (pos_loop.stem_start + 4)):
                        
                        strand = "-" if pos_loop.reverse_compliment else "+"
                        
                        # Add stem-loop feature
                        gff_entry = (f"{record.id}\tsl_finder\tstem_loop\t{pos_loop.stem_start}\t"
                                   f"{pos_loop.stem_end}\t.\t{strand}\t.\tName=stem-loop")
                        features_out.append(gff_entry)
                        
                        # Add nonanucleotide feature
                        nona_entry = (f"{record.id}\tsl_finder\tnonanucleotide\t{pos_loop.motif_start}\t"
                                    f"{pos_loop.motif_start + 8}\t.\t{strand}\t.\tName=nona")
                        features_out.append(nona_entry)
                        
                        # Add to CSV data
                        csv_rows.append({
                            'seqID': record.id,
                            'matched': pos_loop.motif,
                            'motif_start': pos_loop.motif_start,
                            'stem_start': pos_loop.stem_start,
                            'stem_end': pos_loop.stem_end,
                            'score': pos_loop.score,
                            'folded_structure': pos_loop.folded
                        })
                        
                        logger.info(f"Found stem-loop at {pos_loop.stem_start}-{pos_loop.stem_end} in {record.id}")
                
                except Exception as e:
                    logger.error(f"Error processing stem loop in {record.id}: {str(e)}")
                    continue
        
        # Write output files
        out_gff_path = os.path.join(output_dir, out_gff)
        with open(out_gff_path, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            for feat in features_out:
                out_f.write(feat + "\n")
        logger.info(f"Wrote GFF annotations to {out_gff_path}")
        
        if out_csv:
            out_csv_path = os.path.join(output_dir, out_csv)
            with open(out_csv_path, 'w', newline='') as csv_f:
                fieldnames = ['seqID', 'matched', 'motif_start', 'stem_start', 
                            'stem_end', 'score', 'folded_structure']
                writer = csv.DictWriter(csv_f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                for row in csv_rows:
                    writer.writerow(row)
            logger.info(f"Wrote CSV output to {out_csv_path}")
        
        logger.info("SL-Finder module completed successfully.")
        
    except Exception as e:
        logger.error(f"Error in run_sl_finder: {str(e)}")
        raise

def main():
    """Main entry point for the sl_finder module when called via cli.py."""
    parser = argparse.ArgumentParser(description="A module for putative stem-loop annotation")
    parser.add_argument("-i", "--input_fasta", required=True, help="Input FASTA file")
    parser.add_argument("--gff_in", required=True, help="Input GFF/GTF file")
    parser.add_argument("--out_gff", required=True, help="Output GFF filename")
    parser.add_argument("--output", help="Path to the output directory (Default: working directory)", default=".")
    parser.add_argument("--csv_out", help="Output CSV filename", default=None)
    parser.add_argument("--motif", help="Conserved motif (default = nantantan)", default="nantantan")
    parser.add_argument("--family", help="CRESS viral family", 
                       choices=['geminiviridae', 'genomoviridae', 'smacoviridae', 
                               'cycloviridae', 'circoviridae', 'general'])
    parser.add_argument("--idealstemlen", "-s", help="Ideal stem length (default = 11)", 
                       type=int, default=11)
    parser.add_argument("--ideallooplen", "-l", help="Ideal loop length (default = 11)", 
                       type=int, default=11)
    parser.add_argument("--frame", "-f", help="Bases around motif for folding (default = 15)", 
                       type=int, default=15)
    
    args = parser.parse_args()

    # Create output directory
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # Validate FASTA file
    def validate_fasta(filename):
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if any(fasta):
                print("FASTA checked.")
                return filename
            else:
                sys.exit("Error: Input file is not in the FASTA format.\n")

    input_fasta = validate_fasta(args.input_fasta)
    
    try:
        run_sl_finder(
            input_fasta, 
            args.gff_in, 
            args.out_gff, 
            output_dir,
            out_csv=args.csv_out,
            motif=args.motif, 
            family=args.family,
            idealstemlen=args.idealstemlen, 
            ideallooplen=args.ideallooplen,
            frame=args.frame
        )
        return 0
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())


    
