require(ape, quietly = TRUE)       # For tree reading and manipulation
require(dendextend, quietly = TRUE) # For tanglegrams
library(DECIPHER)

tree1 = "/fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree_endonuc/sanitized_sequences.fasta.treefile"
tree2  ="/fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2/tree_helic/sanitized_sequences.fasta.treefile"
output_dir="/fs/project/PAS1117/ricardo/ssDNA_test/test_data/output_2"
metadata_path = "/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/metadata.csv"
fig_name="name.pdf"
height = 15
width = 30

t1 = read.tree(tree1)
t2 = read.tree(tree2)

dend1 <- ReadDendrogram(tree1)
dend2 <- ReadDendrogram(tree2)
metadata <- read.csv(metadata_path) |> dplyr::tibble()

# Extract labels from dendrogram on the left
labels <- dend1 %>% set("labels_to_char") %>% labels 

#Using a metadata table with colours create a vector of colours
labels <- as.data.frame(labels)
labels2 <- merge(labels, metadata, by.x="labels", by.y="protein_description", sort=F)
cols <- as.character(labels2$Colours) 

# Make tanglegram
pdf(file.path(output_dir,fig_name), height = height, width = width)
tanglegram(dend1, dend2,
            highlight_distinct_edges = TRUE,
            highlight_branches_col = TRUE,
            common_subtrees_color_branches = TRUE,
            margin_inner = 20,
            lwd = 2)
dev.off()
