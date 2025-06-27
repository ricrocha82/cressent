#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import argparse
import os
import logging

def detect_input_type(df):
    """
    Automatically detect the input table type based on column names
    Returns: 'prosite' or 'motif_table'
    """
    prosite_columns = ['sequence_ac', 'start', 'stop', 'signature_ac', 'prosite_ann']
    motif_columns = ['seqID', 'motif_name', 'motif_id', 'start', 'end', 'matched']
    
    prosite_match = sum(col in df.columns for col in prosite_columns)
    motif_match = sum(col in df.columns for col in motif_columns)
    
    if prosite_match >= 4:
        return 'prosite'
    elif motif_match >= 5:
        return 'motif_table'
    else:
        raise ValueError("Unable to detect input table type. Expected either Prosite or Motif table format.")

def standardize_dataframe(df, input_type):
    """
    Standardize different input formats to a common structure
    """
    if input_type == 'prosite':
        # Map Prosite columns to standard format
        standardized = df.copy()
        standardized['seq_id'] = df['sequence_ac']
        standardized['motif_type'] = df['prosite_ann']
        standardized['motif_id'] = df['signature_ac']
        standardized['position_start'] = df['start']
        standardized['position_end'] = df['stop']
        standardized['motif_length'] = df['stop'] - df['start'] + 1
        if 'matched' not in standardized.columns:
            standardized['matched'] = 'N/A'
        
    elif input_type == 'motif_table':
        # Map Motif table columns to standard format
        standardized = df.copy()
        standardized['seq_id'] = df['seqID']
        standardized['motif_type'] = df['motif_name']
        standardized['motif_id'] = df['motif_id']
        standardized['position_start'] = df['start']
        standardized['position_end'] = df['end']
        standardized['motif_length'] = df['length']
        standardized['matched'] = df['matched']
    
    return standardized

def create_linear_genome_map(df_std, output_dir, input_type):
    """
    Create a linear genome map showing motif positions along each sequence
    """
    save_path = os.path.join(output_dir, f'genome_map_linear_{input_type}.png')
    
    # Get unique sequences
    sequences = df_std['seq_id'].unique()
    n_seqs = len(sequences)
    
    # Create color map for different motif types
    unique_motifs = df_std['motif_type'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_motifs)))
    color_map = dict(zip(unique_motifs, colors))
    
    # Create figure
    fig, axes = plt.subplots(n_seqs, 1, figsize=(15, 2*n_seqs), sharex=False)
    if n_seqs == 1:
        axes = [axes]
    
    for i, seq_id in enumerate(sequences):
        seq_data = df_std[df_std['seq_id'] == seq_id]
        
        # Get sequence length (approximate from max end position)
        seq_length = seq_data['position_end'].max()
        
        # Plot sequence backbone
        axes[i].plot([0, seq_length], [0, 0], 'k-', linewidth=3, alpha=0.3)
        
        # Plot motifs
        for _, row in seq_data.iterrows():
            start, end = row['position_start'], row['position_end']
            motif_type = row['motif_type']
            
            # Create rectangle for motif
            rect = patches.Rectangle((start, -0.2), end-start, 0.4, 
                                   facecolor=color_map[motif_type], 
                                   edgecolor='black', linewidth=0.5,
                                   alpha=0.8)
            axes[i].add_patch(rect)
            
            # Add motif label
            axes[i].text((start+end)/2, 0.3, row['motif_id'], 
                        ha='center', va='bottom', fontsize=8, rotation=45)
        
        # Formatting
        axes[i].set_ylim(-0.5, 1)
        axes[i].set_xlim(0, seq_length)
        axes[i].set_ylabel(f'{seq_id}', rotation=0, ha='right', va='center')
        axes[i].set_xlabel('Position (aa)')
        axes[i].grid(True, alpha=0.3)
    
    # Create legend at the bottom
    legend_elements = [patches.Patch(facecolor=color_map[motif], 
                                   edgecolor='black', label=motif[:50] + '...' if len(motif) > 50 else motif) 
                      for motif in unique_motifs]
    
    # Adjust layout to make room for legend
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    
    # Add legend at the bottom
    plt.figlegend(handles=legend_elements, loc='lower center', 
                 bbox_to_anchor=(0.5, -0.1), fontsize=9, ncol=3)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    logging.info(f"Linear genome map saved to: {save_path}")

def create_motif_density_plot(df_std, output_dir, input_type):
    """
    Create a density plot showing motif distribution patterns
    """
    save_path = os.path.join(output_dir, f'motif_density_{input_type}.png')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Motif count per sequence using seaborn
    motif_counts = df_std.groupby(['seq_id', 'motif_type']).size().reset_index(name='count')
    sns.barplot(data=motif_counts, x='seq_id', y='count', hue='motif_type', ax=ax1)
    ax1.set_title('Motif Counts per Sequence')
    ax1.set_xlabel('Sequence ID')
    ax1.set_ylabel('Number of Motifs')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend().remove()  # Remove legend from individual plot
    
    # Plot 2: Motif position distribution using seaborn
    # Calculate middle positions
    # df_pos = df_std.copy()
    # df_pos['position'] = (df_pos['position_start'] + df_pos['position_end']) / 2
    
    # sns.histplot(data=df_pos, x='position', hue='motif_type', alpha=0.6, 
    #             bins=20, multiple='stack', ax=ax2)
    # ax2.set_title('Motif Position Distribution')
    # ax2.set_xlabel('Position (aa)')
    # ax2.set_ylabel('Frequency')
    # ax2.legend().remove()  # Remove legend from individual plot

    for motif_type in df_std['motif_type'].unique():
        motif_data = df_std[df_std['motif_type'] == motif_type]
        positions = (motif_data['position_start'] + motif_data['position_end']) / 2
        ax2.hist(positions, alpha=0.6, label=motif_type, bins=20)
    ax2.set_title('Motif Position Distribution')
    ax2.set_xlabel('Position (aa)')
    ax2.set_ylabel('Frequency')
    ax2.legend().remove()  # Remove legend from individual plot
    
    # Create single legend at the bottom
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), 
               fontsize=9, ncol=4, title='Motif Type')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    logging.info(f"Motif density plot saved to: {save_path}")

def create_heatmap_view(df_std, output_dir, input_type):
    """
    Create a heatmap showing motif presence across sequences
    """
    save_path = os.path.join(output_dir, f'motif_heatmap_{input_type}.png')
    
    # Create pivot table: sequences vs motif types
    pivot_data = df_std.groupby(['seq_id', 'motif_type']).size().unstack(fill_value=0)
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(pivot_data, annot=True, cmap='YlOrRd', cbar_kws={'label': 'Motif Count'})
    plt.title(f'Motif Distribution Heatmap ({input_type.title()})')
    plt.xlabel('Motif Type')
    plt.ylabel('Sequence ID')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    logging.info(f"Motif heatmap saved to: {save_path}")

def create_detailed_genome_map(df_std, output_dir, input_type):
    """
    Create a detailed genome map with motif information
    """
    save_path = os.path.join(output_dir, f'detailed_genome_map_{input_type}.png')
    
    sequences = df_std['seq_id'].unique()
    
    # Calculate figure dimensions
    n_seqs = len(sequences)
    fig_height = max(8, n_seqs * 1.5)
    
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    # Color mapping for motif types
    unique_motifs = df_std['motif_type'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_motifs)))
    color_map = dict(zip(unique_motifs, colors))
    
    y_positions = {}
    for i, seq_id in enumerate(sequences):
        y_pos = len(sequences) - i - 1
        y_positions[seq_id] = y_pos
        
        seq_data = df_std[df_std['seq_id'] == seq_id]
        seq_length = seq_data['position_end'].max()
        
        # Draw sequence backbone
        ax.plot([0, seq_length], [y_pos, y_pos], 'k-', linewidth=4, alpha=0.3)
        
        # Draw motifs
        for _, row in seq_data.iterrows():
            start, end = row['position_start'], row['position_end']
            motif_type = row['motif_type']
            
            # Draw motif as colored rectangle
            rect = patches.Rectangle((start, y_pos-0.15), end-start, 0.3,
                                   facecolor=color_map[motif_type],
                                   edgecolor='black', linewidth=0.5,
                                   alpha=0.8)
            ax.add_patch(rect)
            
            # Add position labels for motifs
            if end - start > 2:  # Only label longer motifs
                ax.text((start+end)/2, y_pos+0.4, f"{start}-{end}",
                       ha='center', va='bottom', fontsize=7)
        
        # Add sequence ID
        ax.text(-50, y_pos, seq_id, ha='right', va='center', fontsize=10, weight='bold')
    
    # Formatting
    ax.set_ylim(-0.5, len(sequences) - 0.5)
    ax.set_xlim(-100, df_std['position_end'].max() + 50)
    ax.set_xlabel('Position (amino acids)', fontsize=12)
    ax.set_title(f'Detailed Genome Map - {input_type.title()} Analysis', fontsize=14, weight='bold')
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
    logging.info(f"Detailed genome map saved to: {save_path}")

def create_motif_statistics(df_std, input_type):
    """
    Print summary statistics about the motifs
    """
    print(f"=== MOTIF ANALYSIS SUMMARY ({input_type.upper()}) ===\n")
    
    print(f"Total sequences analyzed: {df_std['seq_id'].nunique()}")
    print(f"Total motifs found: {len(df_std)}")
    print(f"Unique motif types: {df_std['motif_type'].nunique()}\n")
    
    print("Motif type frequency:")
    motif_freq = df_std['motif_type'].value_counts()
    for motif, count in motif_freq.items():
        motif_display = motif[:60] + '...' if len(motif) > 60 else motif
        print(f"  {motif_display}: {count}")
    
    print(f"\nAverage motifs per sequence: {len(df_std) / df_std['seq_id'].nunique():.2f}")
    
    print("\nSequence-specific motif counts:")
    seq_counts = df_std.groupby('seq_id').size()
    for seq, count in seq_counts.items():
        print(f"  {seq}: {count} motifs")
    
    if input_type == 'motif_table' and 'matched' in df_std.columns:
        print(f"\nMotif length distribution:")
        length_stats = df_std['motif_length'].describe()
        print(f"  Mean length: {length_stats['mean']:.1f} aa")
        print(f"  Min length: {length_stats['min']:.0f} aa")
        print(f"  Max length: {length_stats['max']:.0f} aa")

def main():
    """
    Main function to create all visualizations with automatic input detection
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Create genome map visualizations from motif data (Prosite or MEME)")
    parser.add_argument("-f", "--input", required=True, help="Input table file (tab-separated)")
    parser.add_argument("-o", "--output", required=True, help="Output directory to save results")
    parser.add_argument("--format", choices=['prosite', 'motif_table', 'auto'], default='auto', 
                        help="Input format (default: auto-detect)")
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Load the data
    print(f"Loading data from: {args.input}")
    df = pd.read_csv(args.input, sep='\t')
    
    # Detect input type
    if args.format == 'auto':
        input_type = detect_input_type(df)
        print(f"Auto-detected input type: {input_type}")
    else:
        input_type = args.format
        print(f"Using specified input type: {input_type}")
    
    # Standardize dataframe
    df_std = standardize_dataframe(df, input_type)
    
    print(f"Saving results to: {args.output}")
    print("Creating genome map visualizations...")
    
    # Create all visualizations
    create_motif_statistics(df_std, input_type)
    create_linear_genome_map(df_std, args.output, input_type)
    create_detailed_genome_map(df_std, args.output, input_type)
    create_motif_density_plot(df_std, args.output, input_type)
    create_heatmap_view(df_std, args.output, input_type)
    
    print("All visualizations completed!")
    print(f"Files saved in: {args.output}")

if __name__ == "__main__":
    main()