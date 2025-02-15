#!/usr/bin/env Rscript

# Load necessary libraries
library(treeio)
library(tidytree)
library(ggtreeExtra)
library(ggtree)
library(tidyr)
library(dplyr)
# library(treedataverse) # tidytree, treeio, ggtree, and ggtreeExtra
# library(readr)
# devtools::install_github("https://github.com/YuLab-SMU/TDbook.git")
# library(TDbook)

# read the tree
tree <- read.tree("/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree/sub_reps_aligned_trimmed_sequences_sanitized_sequences.fasta.treefile")

# read meta_finals (from alignment and build_tree modules)
meta_data_1 <- read.csv("/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/metadata.csv")
meta_data_2 <- read.delim("/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree/sub_reps_aligned_trimmed_sequences_sanitized_name_table.tsv", sep = "\t")
meta_data_2 = meta_data_2 |> rename(protein_description = Original_Name )
meta_final = meta_data_2 |> left_join(meta_data_1, by = 'protein_description')
meta_final = meta_final |> rename(label = Sanitized_Name)

# joing tree with the meta_final 
trda <- tree %>% 
        left_join(meta_final, by=c("label"))

# layout
layout = "rectangular"
branch.length = "branch.length"
open.angle=0
offset = .14
tip_label = "family"

# if layout == "circular" or "unrooted" we need to change the geom_tiplab to geom_tiplab(aes(angle=angle))
if(layout == "circular" | layout == "unrooted" ){
  p = ggtree(trda,
              mapping = NULL,
              layout = layout,
              open.angle = open.angle,
              mrsd = NULL,
              as.Date = FALSE,
              yscale = "none",
              yscale_mapping = NULL,
              ladderize = TRUE,
              right = FALSE,
              branch.length = branch.length,
              root.position = 0,
              xlim = NULL,
              layout.params = list(),
              hang = 0.1) +
            theme_tree2() +
            geom_tiplab(offset = offset, aes(angle=angle)) +
            geom_tippoint(
                  mapping = aes(
                    shape = tip_label, 
                    color = tip_label
                  )
                ) 
}else{
  p = ggtree(trda,
              mapping = NULL,
              # layout = layout,
              open.angle = open.angle,
              mrsd = NULL,
              as.Date = FALSE,
              yscale = "none",
              yscale_mapping = NULL,
              ladderize = TRUE,
              right = FALSE,
              branch.length = branch.length,
              root.position = 0,
              xlim = NULL,
              layout.params = list(),
              hang = 0.1) +
            theme_tree2() +
            geom_tiplab(offset = offset) +
            geom_tippoint(
                  mapping = aes(
                    shape = tip_label, 
                    color = tip_label
                  )
                ) 
}

# plot tree
# Finally, add tip labels and adjust axis
# p = ggtree(trda,
#               mapping = NULL,
#               layout = layout,
#               open.angle = open.angle,
#               mrsd = NULL,
#               as.Date = FALSE,
#               yscale = "none",
#               yscale_mapping = NULL,
#               ladderize = TRUE,
#               right = FALSE,
#               branch.length = branch.length,
#               root.position = 0,
#               xlim = NULL,
#               layout.params = list(),
#               hang = 0.1) +
#             theme_tree2() +
#             geom_tiplab(offset = offset) +
#             geom_tippoint(
#                   mapping = aes(
#                     shape = tip_label, 
#                     color = tip_label
#                   )
#                 ) 

if(with_alignment == TRUE){
  # phylogenetic tree with aligments
  AA_sequence = "/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree/sub_reps_aligned_trimmed_sequences_sanitized_sequences.fasta"
  p = msaplot(p, fasta=AA_sequence) + theme(legend.position = "none")
}else{
  return(p)
}

ggsave(filename = f, plot = p, width=7, height=7)

# color all the tree
# change the names of the lables to color
new_lables = trda %>%
                select(label, family, isTip)

grp <- split(new_lables$label[new_lables$isTip], new_lables$family[new_lables$isTip])

# plot
p_circ = ggtree(tree, layout = 'circular', branch.length='none')
p_circ = groupOTU(p_circ, grp, 'Family') +
              aes(color=Family) +
              theme(legend.position="right") +
              geom_tiplab(offset = .14)

ggsave(filename = f, plot = p_circ, width=7, height=7)

# Use distance matrix
# here if the argument is use_dist_matrix (--use_dist_matrix)
# Read the file, skipping the first line (which contains "200")
mldist_matrix <- read.table("/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree/sub_reps_aligned_trimmed_sequences_sanitized_sequences.fasta.mldist", skip = 1, header = FALSE, row.names = 1)

# Convert the data frame to a dist matrix
mldist_matrix <- as.dist(mldist_matrix)

# make the tree
tree_color = ape::bionj(mldist_matrix)

# plot the phylogenetic tree
p_circ <- ggtree(tree_color, layout = 'circular', branch.length='none')
p_circ <- groupOTU(p_circ, grp, 'Family') + 
            aes(color=Family) +
            theme(legend.position="right") +
            geom_tiplab(offset = 0)

ggsave(filename = f, plot = p_circ, width=7, height=7)

# plot tree with aligment
# AA_sequence = "/fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output/tree/sub_reps_aligned_trimmed_sequences_sanitized_sequences.fasta"
# p <- ggtree(trda, layout='circular') + 
#     geom_tiplab(offset=4, align=TRUE) + xlim(NA, 12)
# msaplot(p, AA_sequence, window=c(120, 200))

# p <- ggtree(trda) + geom_tiplab(size=3)
# align_p = msaplot(p, AA_sequence, offset=3) + theme(legend.position = "none")

