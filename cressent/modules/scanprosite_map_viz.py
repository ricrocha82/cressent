import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
from collections import defaultdict
import argparse
import os

def load_prosite_data(file_path):
    """Load and prepare Prosite data"""
    df = pd.read_csv(file_path, sep='\t')
    return df

def create_linear_genome_map(df, output_dir):
    """
    Create a linear genome map showing motif positions along each sequence
    """
    save_path = os.path.join(output_dir, 'genome_map_linear.png')
    # Get unique sequences
    sequences = df['sequence_ac'].unique()
    n_seqs = len(sequences)
    
    # Create color map for different motif types
    unique_motifs = df['prosite_ann'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_motifs)))
    color_map = dict(zip(unique_motifs, colors))
    
    # Create figure
    fig, axes = plt.subplots(n_seqs, 1, figsize=(15, 2*n_seqs), sharex=False)
    if n_seqs == 1:
        axes = [axes]
    
    for i, seq_id in enumerate(sequences):
        seq_data = df[df['sequence_ac'] == seq_id]
        
        # Get sequence length (approximate from max stop position)
        seq_length = seq_data['stop'].max()
        
        # Plot sequence backbone
        axes[i].plot([0, seq_length], [0, 0], 'k-', linewidth=3, alpha=0.3)
        
        # Plot motifs
        for _, row in seq_data.iterrows():
            start, stop = row['start'], row['stop']
            motif_type = row['prosite_ann']
            
            # Create rectangle for motif
            rect = patches.Rectangle((start, -0.2), stop-start, 0.4, 
                                   facecolor=color_map[motif_type], 
                                   edgecolor='black', linewidth=0.5,
                                   alpha=0.8)
            axes[i].add_patch(rect)
            
            # Add motif label
            axes[i].text((start+stop)/2, 0.3, row['signature_ac'], 
                        ha='center', va='bottom', fontsize=8, rotation=45)
        
        # Formatting
        axes[i].set_ylim(-0.5, 1)
        axes[i].set_xlim(0, seq_length)
        axes[i].set_ylabel(f'{seq_id}', rotation=0, ha='right', va='center')
        axes[i].set_xlabel('Position (aa)')
        axes[i].grid(True, alpha=0.3)
    
    # Create legend
    legend_elements = [patches.Patch(facecolor=color_map[motif], 
                                   edgecolor='black', label=motif) 
                      for motif in unique_motifs]

    # Adjust layout to make room for legend
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)

    # Add legend at the bottom
    plt.figlegend(handles=legend_elements, loc='lower center', 
                 bbox_to_anchor=(0.5, -0.1), fontsize=9, ncol=3)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def create_motif_density_plot(df, output_dir):
    """
    Create a density plot showing motif distribution patterns
    """
    save_path = os.path.join(output_dir, 'motif_density.png')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Motif count per sequence using seaborn
    motif_counts = df.groupby(['sequence_ac', 'prosite_ann']).size().reset_index(name='count')
    sns.barplot(data=motif_counts, x='sequence_ac', y='count', hue='prosite_ann', ax=ax1)
    ax1.set_title('Motif Counts per Sequence')
    ax1.set_xlabel('Sequence ID')
    ax1.set_ylabel('Number of Motifs')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend().remove()  # Remove legend from individual plot
    
    # Plot 2: Motif position distribution using seaborn
    for motif_type in df['prosite_ann'].unique():
        motif_data = df[df['prosite_ann'] == motif_type]
        positions = (motif_data['start'] + motif_data['stop']) / 2
        ax2.hist(positions, alpha=0.6, label=motif_type, bins=20)
    ax2.set_title('Motif Position Distribution')
    ax2.set_xlabel('Position (aa)')
    ax2.set_ylabel('Frequency')
    ax2.legend().remove()  # Remove legend from individual plot
    
    # Create single legend at the bottom
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1), 
                fontsize=9, ncol=4, title='Motif Type')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def create_heatmap_view(df, output_dir):
    """
    Create a heatmap showing motif presence across sequences
    """
    save_path = os.path.join(output_dir, 'motif_heatmap.png')
    # Create pivot table: sequences vs motif types
    pivot_data = df.groupby(['sequence_ac', 'prosite_ann']).size().unstack(fill_value=0)
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(pivot_data, annot=True, cmap='YlOrRd', cbar_kws={'label': 'Motif Count'})
    plt.title('Motif Distribution Heatmap')
    plt.xlabel('Motif Type')
    plt.ylabel('Sequence ID')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def create_detailed_genome_map(df, output_dir):
    """
    Create a detailed genome map with motif information
    """
    save_path = os.path.join(output_dir, 'detailed_genome_map.png')
    sequences = df['sequence_ac'].unique()
    
    # Calculate figure dimensions
    n_seqs = len(sequences)
    fig_height = max(8, n_seqs * 1.5)
    
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    # Color mapping for motif types
    unique_motifs = df['prosite_ann'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_motifs)))
    color_map = dict(zip(unique_motifs, colors))
    
    y_positions = {}
    for i, seq_id in enumerate(sequences):
        y_pos = len(sequences) - i - 1
        y_positions[seq_id] = y_pos
        
        seq_data = df[df['sequence_ac'] == seq_id]
        seq_length = seq_data['stop'].max()
        
        # Draw sequence backbone
        ax.plot([0, seq_length], [y_pos, y_pos], 'k-', linewidth=4, alpha=0.3)
        
        # Draw motifs
        for _, row in seq_data.iterrows():
            start, stop = row['start'], row['stop']
            motif_type = row['prosite_ann']
            
            # Draw motif as colored rectangle
            rect = patches.Rectangle((start, y_pos-0.15), stop-start, 0.3,
                                   facecolor=color_map[motif_type],
                                   edgecolor='black', linewidth=0.5,
                                   alpha=0.8)
            ax.add_patch(rect)
            
            # Add position labels for important motifs
            if stop - start > 2:  # Only label longer motifs
                ax.text((start+stop)/2, y_pos+0.4, f"{start}-{stop}",
                       ha='center', va='bottom', fontsize=7)
        
        # Add sequence ID
        ax.text(-50, y_pos, seq_id, ha='right', va='center', fontsize=10, weight='bold')
    
    # Formatting
    ax.set_ylim(-0.5, len(sequences) - 0.5)
    ax.set_xlim(-100, df['stop'].max() + 50)
    ax.set_xlabel('Position (amino acids)', fontsize=12)
    ax.set_title('Detailed Genome Map - Prosite Motifs', fontsize=14, weight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Remove y-axis ticks
    ax.set_yticks([])
    
    # Create legend
    legend_elements = [patches.Patch(facecolor=color_map[motif], 
                                   edgecolor='black', label=motif[:40] + '...' if len(motif) > 40 else motif) 
                      for motif in unique_motifs]
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def create_motif_statistics(df):
    """
    Print summary statistics about the motifs
    """
    print("=== MOTIF ANALYSIS SUMMARY ===\n")
    
    print(f"Total sequences analyzed: {df['sequence_ac'].nunique()}")
    print(f"Total motifs found: {len(df)}")
    print(f"Unique motif types: {df['prosite_ann'].nunique()}\n")
    
    print("Motif type frequency:")
    motif_freq = df['prosite_ann'].value_counts()
    for motif, count in motif_freq.items():
        print(f"  {motif}: {count}")
    
    print(f"\nAverage motifs per sequence: {len(df) / df['sequence_ac'].nunique():.2f}")
    
    print("\nSequence-specific motif counts:")
    seq_counts = df.groupby('sequence_ac').size()
    for seq, count in seq_counts.items():
        print(f"  {seq}: {count} motifs")

# Main execution function
def main():
    parser = argparse.ArgumentParser(description="Make a genome map with the Scanprosite results")
    parser.add_argument("-f","--file", help="Input file from scanprosite", required=True)
    parser.add_argument("-o", "--output", default=".", help="Path to the output directory (Default: working directory)")
    
    args = parser.parse_args()
    # Check if the output directory exists.
    if os.path.exists(args.output):
        print(f"Warning: Output directory '{args.output}' already exists. Its contents will be overwritten.")
    else:
        os.makedirs(args.output, exist_ok=True)

    
    # Load the data (replace with your file path)
    df = pd.read_csv(args.file, sep='\t')  # Adjust path as needed
    
    print("Creating genome map visualizations...")
    
    # Create all visualizations
    create_motif_statistics(df)
    create_linear_genome_map(df, args.output)
    create_detailed_genome_map(df, args.output)
    create_motif_density_plot(df, args.output)
    create_heatmap_view(df, args.output)
    
    print("All visualizations completed!")

# Example usage:
if __name__ == "__main__":
    main()