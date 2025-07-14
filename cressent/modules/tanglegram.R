#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(ape)       # For tree reading and manipulation
  library(dendextend) # For tanglegrams
  library(DECIPHER)
})


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
    width = 20,
    lab.cex=1.5
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
lab.cex <- as.numeric(arg_list$lab.cex)

# Load the trees
dend1 <- ReadDendrogram(path1)
dend2 <- ReadDendrogram(path2)

# Extract labels from dendrogram on the left
labels1 <- dend1 %>% set("labels_to_char") %>% labels 
labels2 <- dend2 %>% set("labels_to_char") %>% labels 

# Check for common labels and prune trees if needed
common_labels <- intersect(labels1, labels2)
if (length(common_labels) == 0) {
  stop("Error: The two trees have no common labels.")
} 

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create tanglegram
# Calculate the maximum number of characters in the tree labels
max_label_length <- max(
  max(nchar(labels1)), 
  max(nchar(labels2))
)
# Dynamically set margin_inner based on label length
margin_inner <- max(max_label_length * .5)
# Adjust the Margins
num_labels <- max(length(labels1), length(labels2))

# createa dendrogram list
dl <- dendlist(dend1, dend2)

# # calculate RF score (Robinson-Foulds distance)
RF_score = dist.dendlist(dl, method = c("edgeset"))
print(RF_score)
# save the plot
# Dynamically Adjust the PDF Dimensions
pdf(file = file.path(output_dir, output_name), width = max(width, num_labels / 10), height = max(height, num_labels / 15))
# pdf(file = file.path(output_dir, output_name), width = width, height = height)
dl %>% 
    # dendextend::untangle(method = "step2side") %>%
    tanglegram(main_left = label1,
            main_right = label2,
            main = paste("RF score:", RF_score),
            highlight_distinct_edges = TRUE,
            highlight_branches_col = TRUE,
            common_subtrees_color_branches = TRUE,
            highlight_branches_lwd=TRUE,
            # margin_inner = margin_inner,
            lwd = 2,
            lab.cex=lab.cex)
# dev.copy(pdf, file.path(output_dir, output_name))
dev.off()


cat(sprintf("Tanglegram saved to %s\n", file.path(output_dir, output_name)))

