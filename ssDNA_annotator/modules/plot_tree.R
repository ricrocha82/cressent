#!/usr/bin/env Rscript


# Load required libraries
suppressPackageStartupMessages({
  library(treeio)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ape)
  library(rlang)
})

# --- Parse command-line arguments using commandArgs() ---

args <- commandArgs(trailingOnly = TRUE)

# Set default options in a list
opt <- list(
  tree = NULL,
  outdir = NULL,
  metadata_1 = NULL,
  metadata_2 = NULL,
  alignment = NULL,
  layout = "rectangular",
  branch_length = "branch.length",
  open_angle = 0,
  offset = 0.14,
  tip_label = "family",
  color = FALSE,
  fig_width = 7,
  fig_height = 7,
  plot_tips = TRUE,
  use_dist_matrix = FALSE,
  plot_name = "tree_plot.pdf"
)

# Helper: parse arguments in the form --option=value
for (arg in args) {
  if (grepl("^--", arg)) {
    splitArg <- strsplit(sub("^--", "", arg), "=")[[1]]
    if (length(splitArg) == 2) {
      opt[[splitArg[1]]] <- splitArg[2]
    } else {
      # If flag provided without "=value", set to TRUE
      opt[[splitArg[1]]] <- TRUE
    }
  }
}

# Convert numeric arguments
opt$open_angle <- as.numeric(opt$open_angle)
opt$offset <- as.numeric(opt$offset)
opt$fig_width <- as.numeric(opt$fig_width)
opt$fig_height <- as.numeric(opt$fig_height)

# Convert boolean arguments (if provided as strings)
opt$color <- tolower(as.character(opt$color)) %in% c("true", "1", "yes")
opt$plot_tips <- tolower(as.character(opt$plot_tips)) %in% c("true", "1", "yes")
opt$use_dist_matrix <- tolower(as.character(opt$use_dist_matrix)) %in% c("true", "1", "yes")

# Validate required arguments
if (is.null(opt$tree)) {
  stop("Error: You must provide a tree file using --tree=FILE", call.=FALSE)
}
if (is.null(opt$outdir)) {
  stop("Error: You must provide an output directory using --outdir=DIR", call.=FALSE)
}

# --- Process the tree input ---

if (opt$use_dist_matrix) {
  # Assume the distance matrix file is the tree file name with a .mldist extension
  dist_file <- sub("\\.treefile$", ".mldist", opt$tree)
  message("Reading distance matrix from: ", dist_file)
  mldist_matrix <- read.table(dist_file, skip = 1, header = FALSE, row.names = 1)
  mldist_matrix <- as.dist(mldist_matrix)
  tree <- ape::bionj(mldist_matrix)
} else {
  tree <- read.tree(opt$tree)
}

# --- Merge tree with metadata if provided ---
if (!is.null(opt$metadata_2)) {
  # Both metadata files are provided.
  meta_data_1 <- read.csv(opt$metadata_1)
  meta_data_2 <- read.delim(opt$metadata_2, sep="\t")
  
  meta_data_2 <- meta_data_2 %>% rename(protein_description = Original_Name)
  
  # Merge using protein_description as key, then rename the join column to "label"
  meta_final <- left_join(meta_data_2, meta_data_1, by = "protein_description") %>%
                rename(label = Sanitized_Name)
  
  trda <- tree %>% left_join(meta_final, by = "label")
  
} else if (!is.null(opt$metadata_1)) {
  # Only one metadata file is provided.
  meta_data_1 <- read.csv(opt$metadata_1)
  
  meta_final <- meta_data_1 %>% rename(label = protein_description)
  trda <- tree %>% left_join(meta_final, by = "label")
  
} else {
  trda <- tree
}

# --- Check if the join was successful ---
# Here we assume that if metadata was provided, the 'tip_label' column is expected.
expected_column <- opt$tip_label
cols_to_check = trda %>% select(everything()) |> colnames()

if (!is.null(opt$metadata_1) | !is.null(opt$metadata_2)) {
    if (!(expected_column %in% cols_to_check)) {
    warning(paste("\n\nAfter joining, the expected column", expected_column, 
              "was not found. Please verify that the metadata values match the tree labels.\n\n"))
  }
}

if (!is.null(opt$metadata_1) | !is.null(opt$metadata_2)) {
  missing_matches <- trda %>% select(all_of(expected_column), label) %>% drop_na(all_of(expected_column))
  if(nrow(missing_matches) == 0) {
    warning(paste("\n\nTree tip labels do not have matching metadata in the", 
              expected_column, "column\n\n
              The plot may not be colored adequately.\n\n"))
  }
}


# --- Set plotting parameters ---
layout <- opt$layout
branch_length <- opt$branch_length
open_angle <- opt$open_angle
offset <- opt$offset
tip_label <- opt$tip_label

# --- Create base tree plot ---
if (layout %in% c("circular", "unrooted")) {
  p <- ggtree(trda,
              layout = layout,
              open.angle = open_angle,
              branch.length = branch_length,
              ladderize = TRUE,
              right = FALSE,
              root.position = 0,
              hang = 0.1) +
       theme_tree2()
  if (opt$plot_tips) {
    p <- p + geom_tiplab(offset = offset, aes(angle = angle))
  }
  p <- p + geom_tippoint(mapping = aes(shape = !!sym(tip_label), color = !!sym(tip_label)))
} else {
  p <- ggtree(trda,
              layout = layout,
              open.angle = open_angle,
              branch.length = branch_length,
              ladderize = TRUE,
              right = FALSE,
              root.position = 0,
              hang = 0.1) +
       theme_tree2()
  if (opt$plot_tips) {
    p <- p + geom_tiplab(offset = offset)
  }
  p <- p + geom_tippoint(mapping = aes(shape = !!sym(tip_label), color = !!sym(tip_label)))
}

# --- Add alignment plot if alignment file is provided ---
if (!is.null(opt$alignment)) {
  p <- msaplot(p, fasta = opt$alignment) + theme(legend.position = "none")
}

# --- Optionally color the tree if --color is TRUE and metadata is provided ---
if (opt$color && !is.null(opt$metadata_1) && !is.null(opt$metadata_2)) {
  new_labels <- trda %>% select(label, !!as.name(tip_label), isTip)
  grp <- split(new_labels$label[new_labels$isTip], new_labels[[tip_label]][new_labels$isTip])
  p <- groupOTU(p, grp, tip_label) +
      aes(color = !!sym(tip_label)) +
      theme(legend.position = "right")
}

# --- Save the plot ---
output_file <- file.path(opt$outdir, opt$plot_name)
message("Saving figure to: ", output_file)
ggsave(filename = output_file, plot = p, width = opt$fig_width, height = opt$fig_height)

# ./ssDNA_annotator/modules/plot_tree.R --tree=mytree.treefile --outdir=output_dir \
#                   --metadata_1=meta.csv --metadata_2=metadata_2.tsv \
#                   --alignment=alignment.fasta \
#                   --layout=circular --branch_length=branch.length \
#                   --open_angle=90 --offset=0.2 --tip_label=family 
#                   --color=TRUE --fig_width=8 --fig_height=8 --plot_tips=TRUE --use_dist_matrix=FALSE \
#                   --plot_name=my_custom_tree.pdf