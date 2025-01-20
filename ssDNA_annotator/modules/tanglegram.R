#!/usr/bin/env Rscript

# Load necessary libraries
require(ape, quietly = TRUE)       # For tree reading and manipulation
require(dendextend, quietly = TRUE) # For tanglegrams

# Define CLI arguments
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Define default values for arguments
arg_list <- list(
    tree1 = NULL,
    tree2 = NULL,
    label1 = "Tree 1",
    label2 = "Tree 2",
    output = "./",
    name_tanglegram = "tanglegram.pdf",
    height = 11,
    width = 20
)

# Assign values from command-line arguments
for (arg in args) {
  key_value <- strsplit(arg, "=")[[1]]
  if (length(key_value) == 2) {
    arg_list[[key_value[1]]] <- key_value[2]
  }
}

# Check for required arguments
if (is.null(arg_list$tree1) || is.null(arg_list$tree2)) {
  stop("Error: --tree1 and --tree2 are required arguments.")
}


# Assign arguments to variables
path1 <- arg_list$tree1
path2 <- arg_list$tree2
label1 <- arg_list$label1
label2 <- arg_list$label2
output_dir <- arg_list$output
output_name <- arg_list$name_tanglegram
height <- as.numeric(arg_list$height)
width <- as.numeric(arg_list$width)

# Load the trees
tree1 <- ape::read.tree(path1)
tree2 <- ape::read.tree(path2)

# Check for common labels and prune trees if needed
common_labels <- intersect(tree1$tip.label, tree2$tip.label)
if (length(common_labels) == 0) {
  stop("Error: The two trees have no common labels.")
} else {
    tree1_labels <- tree1$tip.label
    tree2_labels <- tree2$tip.label
}

# Generate distance matrices
tree1_dists <- cophenetic(tree1)
tree2_dists <- cophenetic(tree2)

# Generate dendrograms
tree1_hc <- hclust(as.dist(tree1_dists), method = "complete")
tree1_dend <- as.dendrogram(tree1_hc)
tree2_hc <- hclust(as.dist(tree2_dists), method = "complete")
tree2_dend <- as.dendrogram(tree2_hc)

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create tanglegram
# Calculate the maximum number of characters in the tree labels
max_label_length <- max(
  max(nchar(tree1$tip.label)), 
  max(nchar(tree2$tip.label))
)
# Dynamically set margin_inner based on label length
margin_inner <- max(max_label_length * 0.5)
# Adjust the Margins
num_labels <- max(length(tree1$tip.label), length(tree2$tip.label))

# createa dendrogram list
dl <- dendlist(tree1_dend, tree2_dend)

# save the plot
# Dynamically Adjust the PDF Dimensions
pdf(file = file.path(output_dir, output_name), width = max(width, num_labels / 10), height = max(height, num_labels / 15))
# pdf(file = file.path(output_dir, output_name), width = width, height = height)
dl %>% 
    # dendextend::untangle(method = "step2side") %>%
    tanglegram(main_left = label1,
            main_right = label2,
            sort = TRUE,
            highlight_distinct_edges = TRUE,
            highlight_branches_col = TRUE,
            common_subtrees_color_branches = TRUE,
            margin_inner = margin_inner,
            lwd = 2)
# dev.copy(pdf, file.path(output_dir, output_name))
dev.off()


cat(sprintf("Tanglegram saved to %s\n", file.path(output_dir, output_name)))

# Rscript /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/tanglegram.R \
#     tree1=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/tree_helic/sanitized_sequences.fasta.treefile \
#     tree2=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/tree_endonuc/sanitized_sequences.fasta.treefile \
#     label1="Helicase" \
#     label2="Endonuclease" \
#     output=/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_2/ \
#     name_tanglegram="my_tanglegram.pdf" \
#     height=11 \
#     width=20