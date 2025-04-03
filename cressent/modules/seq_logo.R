suppressPackageStartupMessages({
  require(ggplot2, quietly = TRUE)
  require(ggseqlogo, quietly = TRUE)
  require(DECIPHER, quietly = TRUE)
  require(Biostrings, quietly = TRUE)
  require(dplyr, quietly = TRUE)
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
    group_label=NULL
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

cat("Output directory:", output_dir, "\n")
cat("Output name:", output_name, "\n")

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

# Read input sequences
if (!is.null(fasta_file) && file.exists(fasta_file)) {
  seqs <- readBStringSet(fasta_file)
  seq_type <- detect_sequence_type(as.character(seqs))
} else if (!is.null(seq_df) && file.exists(seq_df)) {
  df_loc <- read.delim(seq_df, sep = "\t")
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

# Align sequences if necessary
if (length(unique(width(seqs))) > 1) {
  seqs <- AlignSeqs(seqs)
}

# Convert to character for ggseqlogo
motifs <- as.character(seqs)

# Generate sequence logo
p_mot <- ggplot() + geom_logo(motifs) + theme_logo() + labs(title = plot_title)

# Run this if the user wants to split the sequence logo by a determined group
if (split && !is.null(metadata) && !is.null(ncol) && !is.null(group_label)) {
  meta_data <- read.csv(metadata)
  merged_df <- df_loc %>% left_join(meta_data, by = join_by(seqID == protein_description))

  # Split data based on the specified group
  grp <- split(merged_df$matched, merged_df[[group_label]])
  
  # Generate a grouped sequence logo
  p_mot <- ggseqlogo(grp, ncol = ncol) + theme_logo()
}

# Save plot
path_out <- file.path(output_dir, output_name)
ggsave(p_mot, filename = path_out, width = width, height = height)


# Rscript /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/seq_logo.R \
#     seq_df=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/motif/pattern_positions.txt \
#     output_dir=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/motif \
#     output_name="my_seq_log_test.pdf" \
#     plot_title=test \
#     height=10 \
#     width=5 \
#     split=TRUE \
#     metadata=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/metadata.csv






