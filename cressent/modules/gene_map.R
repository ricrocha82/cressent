#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
library(ggplot2)
library(gggenes)
})


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage information
usage <- function() {
  cat("Usage: Rscript gene_map.R [input_file] [output_file] [height] [width] [plot_title]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  input_file    Path to the motif table CSV (tab-delimited)\n")
  cat("  output_file   Path for the output PDF file (default: gene_motif.pdf)\n")
  cat("  height        Height of the output plot in inches (default: 10)\n")
  cat("  width         Width of the output plot in inches (default: 10)\n")
  cat("  plot_title    Title for the plot (optional)\n")
  cat("\n")
  cat("Example: Rscript gene_map.R motif_table.csv gene_motif.pdf 8 12 \"Motif Distribution\"\n")
  quit(status = 1)
}

# Check if help is requested
if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  usage()
}

# Set default values
input_file <- NULL
output_file <- "gene_map_motif"
height <- 10
width <- 10
plot_title <- NULL

# Parse arguments
if (length(args) >= 1) input_file <- args[1]
if (length(args) >= 2) output_file <- args[2]
if (length(args) >= 3) height <- as.numeric(args[3])
if (length(args) >= 4) width <- as.numeric(args[4])
if (length(args) >= 5) plot_title <- args[5]

# Check if input file is provided
if (is.null(input_file)) {
  cat("Error: Input file must be specified.\n\n")
  usage()
}


# Check if the output directory exists, create it if it doesn't
output_dir <- dirname(output_file)
if (!dir.exists(output_dir) && output_dir != ".") {
  cat(paste("Creating output directory:", output_dir, "\n"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    cat(paste("Error: Failed to create output directory:", output_dir, "\n"))
    quit(status = 1)
  }
}


# Read and process the data
cat(paste("Reading data from:", input_file, "\n"))
df_genes <- try(read.delim(input_file, sep = "\t"))
if (inherits(df_genes, "try-error")) {
  cat("Error: Failed to read input file. Check if it's a valid tab-delimited file.\n")
  quit(status = 1)
}

# Count rows for legend
n_rows <- length(unique(df_genes$matched))
cat(paste("Found", n_rows, "unique motifs\n"))

# Create the plot
cat("Generating plot...\n")
plot_gene <- ggplot(df_genes, aes(xmin = start, xmax = end, y = seqID, fill = matched)) +
                geom_gene_arrow() +
                facet_wrap(~ seqID, scales = "free", ncol = 1) +
                scale_fill_brewer(palette = "Set3") +
                guides(fill=guide_legend(nrow=n_rows, byrow=TRUE)) +
                labs(y = "Sequence name", fill = "Motifs") +
                theme_genes() +
                theme(legend.position = "bottom")

# Add title if provided
if (!is.null(plot_title)) {
  cat(paste("Adding plot title:", plot_title, "\n"))
  plot_gene <- plot_gene + ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

# Save the plot
cat(paste("Saving plot to:", output_file, "\n"))
ggsave(plot_gene, filename = file.path(output_file), height = height, width = width)
cat("Done!\n")

# Rscript /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/gene_map.R \
#                     /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc/motif_table.csv \
#                     /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/motif_disc/gene_map_motif.pdf \
#                     10 10 \
#                     "Motif Distribution in Genes"