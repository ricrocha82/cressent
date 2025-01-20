require(ggplot2, quietly = TRUE)
require(ggseqlogo, quietly = TRUE)
require(DECIPHER, quietly = TRUE)
require(Biostrings, quietly = TRUE)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

arg_list <- list(
    seq_df = NULL,
    plot_title = 'sequence_logo',
    output_dir = "./",
    output_name = "sequence_logo.pdf",
    width = 10, 
    height = 5

)

# Assign values from command-line arguments
for (arg in args) {
  key_value <- strsplit(arg, "=")[[1]]
  if (length(key_value) == 2) {
    arg_list[[key_value[1]]] <- key_value[2]
  }
}

# Check for required arguments
if (is.null(arg_list$seq_df)) {
  stop("Error: --motif_df is a required argument.")
}

# Assign arguments to variables
path1 <- arg_list$seq_df
plot_title <- arg_list$plot_title
output_dir <- arg_list$output_dir
output_name <- arg_list$output_name
width <- as.numeric(arg_list$width)
height <- as.numeric(arg_list$height)

cat("Output directory:", output_dir, "\n")
cat("Output name:", output_name, "\n")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# get the args
df_loc = read.delim(path1, sep = "\t")
motifs = df_loc$matched
# check if all the sequences have the same length
# if not, performe the aligment
if(length(unique(motifs)) == 1){
  return(motifs)
} else {
  seqs = DNAStringSet(motifs)
  # Perform multiple sequence alignment
  aligned <- AlignSeqs(seqs)
  motifs = as.character(aligned)
}

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