import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def calculate_gc_content(sequence, window_size=30, step_size=5):
    """
    Calculates %G+C content using a sliding window approach.
    
    Parameters:
        sequence (str): The DNA sequence.
        window_size (int): Size of the sliding window.
        step_size (int): Step size for the sliding window.

    Returns:
        gc_values (list): List of %G+C content for each window.
    """
    gc_values = []
    for i in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[i:i + window_size]
        if "N" in window or "-" in window:  # Mark as ambiguous if contains gaps or N's
            gc_values.append(np.nan)
        else:
            gc_count = sum(1 for base in window if base in "GC")
            gc_values.append((gc_count / window_size) * 100)
    return gc_values

def pad_sequences(matrix, fill_value=np.nan):
    """
    Pads sequences in a list of lists to have the same length.
    
    Parameters:
        matrix (list of lists): List containing GC content values.
        fill_value (float): Value to use for padding (default: NaN).
    
    Returns:
        np.array: 2D NumPy array with equal-length rows.
    """
    max_length = max(len(row) for row in matrix)
    return np.array([row + [fill_value] * (max_length - len(row)) for row in matrix])

def plot_gc_heatmap(input_fasta, output_dir=".", window_size=30, step_size=5,
                    xticklabels=10, fig_width=12, fig_height=6, output_name="gc_heatmap.png"):
    """
    Reads a FASTA file and generates a heatmap of %G+C content.
    
    Parameters:
        input_fasta (str): Path to the FASTA file.
        output_dir (str): Directory to save the output file.
        window_size (int): Size of the sliding window.
        step_size (int): Step size for the sliding window.
        xticklabels (int): Interval for x-axis tick labels.
        fig_width (int): Width of the figure.
        fig_height (int): Height of the figure.
        output_name (str): Name of the output image file.
    """
    sequences = []
    seen_names = {}  # Dictionary to track duplicate names
    gc_matrix = []
    avg_gc_content = []

    # Read sequences from FASTA
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_name = record.id

        # Handle duplicate names by appending description
        if seq_name in seen_names:
            seq_name = f"{record.id}_{record.description.replace(' ', '_')}"
        seen_names[seq_name] = True  # Mark name as seen

        seq = str(record.seq).upper()
        gc_values = calculate_gc_content(seq, window_size, step_size)
        
        sequences.append(seq_name)
        gc_matrix.append(gc_values)
        avg_gc_content.append(np.nanmean(gc_values))  # Ignore NaN values for average GC calculation

    # Pad sequences to ensure uniform matrix size
    gc_matrix = pad_sequences(gc_matrix)

    # Define rainbow color map with green as midpoint
    cmap = sns.color_palette("rainbow", as_cmap=True)

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(gc_matrix, cmap=cmap, linewidths=0.5, linecolor='gray',
                xticklabels=xticklabels, 
                yticklabels=sequences, mask=np.isnan(gc_matrix), cbar=True, ax=ax,
                cbar_kws=dict(use_gridspec=False, location="bottom"))  # Legend at the bottom

    # Customize plot
    ax.set_xlabel("Genomic Window")
    ax.set_ylabel("Sequences")

    # Add right-side GC content summary
    for i, avg_gc in enumerate(avg_gc_content):
        ax.text(gc_matrix.shape[1] + 2, i + 0.5, f"{avg_gc:.2f}%", va="center", ha="left", fontsize=10)
        
    # Add labels for right-side table
    plt.text(gc_matrix.shape[1] + 1, -0.5, "Avg. %G+C", fontsize=12, fontweight="bold")

    # Save plot
    output_path = os.path.join(output_dir, output_name)
    plt.savefig(output_path, bbox_inches="tight")
    print(f"GC Heatmap saved to: {output_path}")

    # Close plot to prevent display issues in CLI
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate a GC content heatmap from a FASTA file.")

    parser.add_argument("-i","--input_fasta", help="Input FASTA file containing sequences.")
    parser.add_argument("-o","--output", default=".", help="Path to the output directory (Default: working directory)")
    parser.add_argument("--window_size", type=int, default=30, help="Sliding window size for GC content calculation (default: 30).")
    parser.add_argument("--step_size", type=int, default=5, help="Step size for GC calculation (default: 5).")
    parser.add_argument("--xticklabels", type=int, default=10, help="Interval for x-axis tick labels (default: None).")
    parser.add_argument("--fig_width", type=int, default=12, help="Figure width in inches (default: 12).")
    parser.add_argument("--fig_height", type=int, default=6, help="Figure height in inches (default: 6).")
    parser.add_argument("--output_name", default="gc_heatmap.pdf", help="Name of output image file with extension (default: gc_heatmap.png).")

    args = parser.parse_args()

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
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
    input_fasta = validate_fasta(args.input_fasta)

    # run the script
    plot_gc_heatmap(input_fasta, output_dir, args.window_size, args.step_size,
                    args.xticklabels, args.fig_width, args.fig_height, args.output_name)

if __name__ == "__main__":
    main()

