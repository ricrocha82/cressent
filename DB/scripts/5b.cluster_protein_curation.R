library(dplyr)
library(data.table)
# library(Biostrings) # For working with FASTA files


# Set directory and file paths
dir = '/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db'
family_name = read.delim('/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/familycopy.txt', sep = '\t', header = FALSE)

###################################################3
######### mannually curate the data ################
#####################################################

# create folders to save the annotated and unannotated proteins for each family
dir_path_ann = paste0(dir, "/annotated/")
dir_path_unann = paste0(dir, "/unannotated/")

# Function to create directory if it doesn't exist
create_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created successfully: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}
# Create directories
create_dir(dir_path_ann)
create_dir(dir_path_unann)


# Function to write FASTA files for each cluster
write_fasta_by_cluster <- function(df, family, output_dir) {
    # Ensure the output directory exists
    create_dir(output_dir)
    
    # Group by cluster and write each group's data to a FASTA file
    df %>%
        group_by(cluster) %>%
        group_walk(~ writeLines(
            paste0(">", .x$protein_id, " ", .x$protein_description, "\n", .x$protein_sequence),
            file.path(output_dir, paste0(family,"_",.y$cluster, ".fa"))
            ))
}


process_family <- function(family, dir, dir_path_ann, dir_path_unann) {

  print(family)

  #############################3
  #### Cap proteins ############
  ##############################

  df_longer = fread(paste0(dir,"/mcl/",family,"_protein_description.csv"))

  # filter annotated proteins
  df_caps =  df_longer[
    grepl("cap|Cap|capsid|Capsid|VP", protein_description) &
    !grepl("non-capsid|non-Capsid", protein_description)
  ]

  # get cluster with >=10 annotated Caps or VP proteins
  clusters_to_keep_cap <- df_caps[, .N, by = cluster][N >= 8, cluster]
  # clusters_to_keep_cap = df_caps %>% 
  #                     group_by(cluster) %>%
  #                     count() %>%
  #                     arrange(desc(n)) %>% 
  #                     filter(n >= 10) %>%
  #                     pull(cluster)

  if(length(clusters_to_keep_cap) > 0){
    df_caps = df_caps[cluster %in% clusters_to_keep_cap]
    df_caps_to_save = df_caps %>% select(cluster,protein_id,protein_description)
    create_dir(paste0(dir_path_ann,"caps/"))
    fwrite(df_caps_to_save, file.path(dir_path_ann, "caps", paste0(family, "_caps.csv")))
    # write_csv(df_caps_to_save, file = paste0(dir_path_ann,"caps/" ,family,"_caps.csv"))

    # write fasta
    write_fasta_by_cluster(df_caps, family, output_dir = paste0(dir_path_ann,"caps/"))

  } else {
      clusters_to_keep_cap <- df_caps[, .N, by = cluster][N >= 1, cluster]
      df_caps = df_caps[cluster %in% clusters_to_keep_cap]
      df_caps_to_save = df_caps %>% select(cluster,protein_id,protein_description)
      create_dir(paste0(dir_path_ann,"caps/"))
      fwrite(df_caps_to_save, file.path(dir_path_ann, "caps", paste0(family, "_caps.csv")))
      print(paste("the dataset has less than 8 Caps per cluster:", family))
  } # else {
  #     print(paste("dataset empty for Cap in", family))
  # }

  #############################3
  #### Rep proteins ############
  ##############################
  df_rep <- df_longer[
    grepl("rep|Rep|replication|Replication|replicase|Replicase", protein_description) &
    !grepl("non-rep|non-Rep|non-replication|non-Replication|non-replicase|non-Replicase", protein_description)
  ]
  
   # get cluster with >=8 annotated Reps
  clusters_to_keep_rep <- df_rep[, .N, by = cluster][N >= 1, cluster]
  # clusters_to_keep_rep = df_rep %>% 
                          # group_by(cluster) %>%
                          # count() %>%
                          # arrange(desc(n)) %>%
                          # filter(n >= 8) %>%
                          # pull(cluster)

  if(length(clusters_to_keep_rep) > 0){
      df_rep = df_rep[cluster %in% clusters_to_keep_rep]

      df_rep_to_save = df_rep %>% select(cluster,protein_id,protein_description)
      create_dir(paste0(dir_path_ann,"reps/"))
      fwrite(df_rep_to_save, file.path(dir_path_ann, "reps/", paste0(family, "_reps.csv")))
      # write_csv(df_rep_to_save, file = paste0(dir_path_ann, "reps/",family,"_reps.csv") )

      # write fasta
      write_fasta_by_cluster(df_rep, family, output_dir = paste0(dir_path_ann,"reps/"))

      } else {
        clusters_to_keep_rep <- df_rep[, .N, by = cluster][N >= 1, cluster]
        df_rep = df_rep[cluster %in% clusters_to_keep_rep]
        df_rep_to_save = df_rep %>% select(cluster,protein_id,protein_description)
        create_dir(paste0(dir_path_ann,"reps/"))
        fwrite(df_rep_to_save, file.path(dir_path_ann, "reps/", paste0(family, "_reps.csv")))
          # print(paste("dataset empty for Rep in", family))
  }



  # ###################################
  # #### unannotated ORFs ############
  # ###################################
  df_orf = df_orf <- df_longer[
    !(cluster %in% clusters_to_keep_cap | cluster %in% clusters_to_keep_rep)
  ]

  # get cluster with >=10 unannotated ORF
  clusters_to_keep_orf <- df_orf[, .N, by = cluster][N >= 10, cluster]
  # clusters_to_keep_orf = df_orf %>% 
  #                         group_by(cluster) %>%
  #                         count() %>%
  #                         arrange(desc(n)) %>%
  #                         filter(n >= 10) %>%
  #                         pull(cluster)

  if (length(clusters_to_keep_orf) > 0) {
      df_orf = df_orf[cluster %in% clusters_to_keep_orf]
      create_dir(paste0(dir_path_unann, family,"/"))
      df_orf_to_save = df_orf %>% select(cluster,protein_id,protein_description)
      fwrite(df_orf_to_save, file.path(dir_path_unann, family, paste0(family, "_orf.csv")))
      # write_csv(df_orf_to_save, file = paste0(dir_path_unann, family,"/",family,"_orf.csv") )

      # write fasta
      write_fasta_by_cluster(df_orf, family, output_dir = paste0(dir_path_unann, family))

  } else {
      print(paste("dataset empty for ORF in", family))
  }

}

# Main loop
for (family in family_name$V1) {
  process_family(family, dir, dir_path_ann, dir_path_unann)
}