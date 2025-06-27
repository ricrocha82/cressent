suppressPackageStartupMessages({
  require(ggplot2, quietly = TRUE)
  require(ggseqlogo, quietly = TRUE)
  require(DECIPHER, quietly = TRUE)
  require(Biostrings, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(patchwork, quietly = TRUE)  # For combining plots
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

arg_list <- list(
    fasta_file = NULL,
    seq_df = NULL,
    plot_title = 'sequence_logo',
    output_dir = "./",
    output_name = "sequence_logo.pdf",
    width = 10, 
    height = 5,
    split = FALSE,
    metadata = NULL,
    ncol = NULL,
    group_label = NULL,
    positions_per_row = 50,  # New parameter for positions per row
    max_positions_single_row = 100,  # New parameter for automatic row splitting
    method = 'prob'  # New parameter for ggseqlogo method
)

# Assign values from command-line arguments
for (arg in args) {
  key_value <- strsplit(arg, "=")[[1]]
  if (length(key_value) == 2) {
    arg_list[[key_value[1]]] <- key_value[2]
  }
}

# Assign arguments to variables
fasta_file <- arg_list$fasta_file
seq_df <- arg_list$seq_df
plot_title <- arg_list$plot_title
output_dir <- arg_list$output_dir
output_name <- arg_list$output_name
width <- as.numeric(arg_list$width)
height <- as.numeric(arg_list$height)
split <- as.logical(arg_list$split)
metadata <- arg_list$metadata
ncol <- as.numeric(arg_list$ncol)
group_label <- arg_list$group_label
positions_per_row <- as.numeric(arg_list$positions_per_row)
max_positions_single_row <- as.numeric(arg_list$max_positions_single_row)
method <- arg_list$method

cat("Output directory:", output_dir, "\n")
cat("Output name:", output_name, "\n")
cat("Positions per row:", positions_per_row, "\n")
cat("Max positions for single row:", max_positions_single_row, "\n")
cat("ggseqlogo method:", method, "\n")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to detect sequence type (Protein vs Nucleotide)
detect_sequence_type <- function(seqs) {
  nucleotide_chars <- c("A", "T", "C", "G", "U", "N", "-")
  seq_chars <- unique(strsplit(paste0(seqs, collapse = ""), "")[[1]])
  
  # If sequence contains non-nucleotide characters (e.g., F, L, M, Y), assume it's a protein
  if (any(!seq_chars %in% nucleotide_chars)) {
    return("protein")
  } else {
    return("nucleotide")
  }
}

# Function to split sequences into chunks for multi-row display
split_sequences_into_rows <- function(motifs, positions_per_row) {
  seq_length <- nchar(motifs[1])
  
  if (seq_length <= positions_per_row) {
    return(list(motifs))
  }
  
  # Calculate number of rows needed
  num_rows <- ceiling(seq_length / positions_per_row)
  
  # Split sequences into chunks
  seq_chunks <- list()
  for (i in 1:num_rows) {
    start_pos <- (i - 1) * positions_per_row + 1
    end_pos <- min(i * positions_per_row, seq_length)
    
    chunk_motifs <- substr(motifs, start_pos, end_pos)
    seq_chunks[[i]] <- chunk_motifs
  }
  
  return(seq_chunks)
}

# Function to create multi-row sequence logo
create_multirow_logo <- function(motifs, plot_title, positions_per_row, max_positions_single_row, method) {
  seq_length <- nchar(motifs[1])
  
  # If sequences are short enough, create single plot
  if (seq_length <= max_positions_single_row) {
    return(ggplot() + geom_logo(motifs, method = method) + theme_logo() + 
           labs(title = plot_title) +
           theme(plot.title = element_text(hjust = 0.5)))
  }
  
  # Split sequences into rows
  seq_chunks <- split_sequences_into_rows(motifs, positions_per_row)
  
  # Create individual plots for each row
  plot_list <- list()
  for (i in 1:length(seq_chunks)) {
    start_pos <- (i - 1) * positions_per_row + 1
    end_pos <- min(i * positions_per_row, seq_length)
    chunk_length <- nchar(seq_chunks[[i]][1])
    
    row_title <- if (length(seq_chunks) > 1) {
      paste0("Positions ", start_pos, "-", end_pos)
    } else {
      ""
    }
    
    # Remove legend from all plots except the last one
    show_legend <- (i == length(seq_chunks))
    
    # Calculate breaks for x-axis (every 5 positions, but ensure they fit within chunk)
    x_breaks <- seq(1, chunk_length, by = 5)
    if (max(x_breaks) < chunk_length) {
      x_breaks <- c(x_breaks, chunk_length)
    }
    
    # Calculate corresponding actual position labels
    x_labels <- seq(start_pos, end_pos, by = 5)
    if (length(x_labels) < length(x_breaks)) {
      x_labels <- c(x_labels, end_pos)
    }
    
    # Ensure breaks and labels have the same length
    min_length <- min(length(x_breaks), length(x_labels))
    x_breaks <- x_breaks[1:min_length]
    x_labels <- x_labels[1:min_length]
    
    # Create the plot with correct x-axis positions
    p <- ggplot() + 
      geom_logo(seq_chunks[[i]], method = method) + 
      theme_logo() +
      labs(title = row_title) +
      # Set x-axis breaks and labels to show continuous positions
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels,
        limits = c(0.5, chunk_length + 0.5)
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = if (show_legend) "bottom" else "none",
        plot.margin = margin(t = 5, r = 5, b = if (show_legend) 10 else 2, l = 5, unit = "pt"),
        axis.text.x = element_text(size = 8)
      )
    
    plot_list[[i]] <- p
  }
  
  # Combine plots vertically using patchwork with minimal spacing
  combined_plot <- wrap_plots(plot_list, ncol = 1) &
    theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"))
  
  # Add main title with reduced spacing
  final_plot <- combined_plot + 
    plot_annotation(
      title = plot_title,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
      )
    )
  
  return(final_plot)
}

# Read input sequences
if (!is.null(fasta_file) && file.exists(fasta_file)) {
  seqs <- readBStringSet(fasta_file)
  seq_type <- detect_sequence_type(as.character(seqs))
} else if (!is.null(seq_df) && file.exists(seq_df)) {
  df_loc <- read.delim(seq_df, sep = "\t")
  # Convert all lowercase letters to uppercase in the matched column
  df_loc$matched <- toupper(df_loc$matched)
  seq_type <- detect_sequence_type(df_loc$matched)
  
  if (seq_type == "protein") {
    seqs <- AAStringSet(df_loc$matched)
  } else {
    seqs <- DNAStringSet(df_loc$matched)
  }
} else {
  stop("Error: Either --fasta_file or --seq_df must be provided.")
}

cat("Detected sequence type:", seq_type, "\n")
cat("Number of sequences:", length(seqs), "\n")

# Align sequences if necessary
if (length(unique(width(seqs))) > 1) {
  cat("Aligning sequences...\n")
  seqs <- AlignSeqs(seqs)
}

# Convert to character for ggseqlogo
motifs <- as.character(seqs)
seq_length <- nchar(motifs[1])
cat("Sequence length:", seq_length, "\n")

# Generate sequence logo
if (split && !is.null(metadata) && !is.null(ncol) && !is.null(group_label)) {
  # Handle grouped/split sequence logos
  meta_data <- read.csv(metadata)
  merged_df <- df_loc %>% left_join(meta_data, by = join_by(seqID == protein_description))
  
  # Split data based on the specified group
  grp <- split(merged_df$matched, merged_df[[group_label]])
  
  # For grouped data, we'll create separate multi-row plots for each group
  if (seq_length > max_positions_single_row) {
    # Create multi-row plots for each group
    group_plots <- list()
    for (group_name in names(grp)) {
      group_motifs <- grp[[group_name]]
      group_plot <- create_multirow_logo(group_motifs, group_name, positions_per_row, max_positions_single_row, method)
      group_plots[[group_name]] <- group_plot
    }
    
    # Combine group plots
    p_mot <- wrap_plots(group_plots, ncol = ncol)
    p_mot <- p_mot + plot_annotation(title = plot_title)
  } else {
    # Use original grouped approach for shorter sequences
    p_mot <- ggseqlogo(grp, method = method, ncol = ncol) + theme_logo() + labs(title = plot_title)
  }
} else {
  # Handle single sequence logo (potentially multi-row)
  p_mot <- create_multirow_logo(motifs, plot_title, positions_per_row, max_positions_single_row, method)
}

# Adjust height based on number of rows
if (seq_length > max_positions_single_row) {
  num_rows <- ceiling(seq_length / positions_per_row)
  # Adjust height proportionally to number of rows with tighter spacing
  height <- height * num_rows * 0.6  # Reduced factor for closer plots
  cat("Adjusted height for", num_rows, "rows:", height, "\n")
}

# Save plot
path_out <- file.path(output_dir, output_name)
ggsave(p_mot, filename = path_out, width = width, height = height, limitsize = FALSE)

cat("Sequence logo saved to:", path_out, "\n")

# Example usage:
# Rscript seq_logo.R \
#     seq_df=pattern_positions.txt \
#     output_dir=output \
#     output_name="multi_row_logo.pdf" \
#     plot_title="Multi-row Sequence Logo" \
#     positions_per_row=50 \
#     max_positions_single_row=100 \
#     method=bits \
#     height=8 \
#     width=12