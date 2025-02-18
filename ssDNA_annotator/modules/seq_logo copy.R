require(ggplot2, quietly = TRUE)
require(ggseqlogo, quietly = TRUE)
require(DECIPHER, quietly = TRUE)
require(Biostrings, quietly = TRUE)

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
    group_label = NULL

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
group_label <- arg_list$group_label

cat("Output directory:", output_dir, "\n")
cat("Output name:", output_name, "\n")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read input sequences
if (!is.null(fasta_file) && file.exists(fasta_file)) {
  seqs <- readBStringSet(fasta_file)
} else if (!is.null(seq_df) && file.exists(seq_df)) {
  df_loc <- read.delim(seq_df, sep = "\t")
  seqs <- DNAStringSet(df_loc$matched)
} else {
  stop("Error: Either --fasta_file or --seq_df must be provided.")
}

# Align sequences if necessary
if (length(unique(width(seqs))) > 1) {
  seqs <- AlignSeqs(seqs)
}

# Convert to character for ggseqlogo
motifs <- as.character(seqs)

# Generate sequence logo
p_mot = ggplot() + geom_logo(motifs) + theme_logo() + labs(title = plot_title)

# save plot
path_out = file.path(output_dir, output_name)
ggsave(p_mot , filename = path_out, width = width , height = height)


# Rscript /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/seq_logo.R \
#     seq_df=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/motif_positions.txt \
#     output_dir=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/ \
#     output_name="my_seq_log_test.pdf" \
#     plot_title=test \
#     height=10 \
#     width=5




seqs_dna
