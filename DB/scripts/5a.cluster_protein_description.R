library(tidyverse)
library(data.table)
library(Biostrings) # For working with FASTA files

# Set directory and file paths
dir = '/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db'
family_name = read_delim('/fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/family.txt', delim = '\t', col_names = FALSE)


# Read family information
for(i in seq_along(family_name$X1)){
    family = as.data.frame(family_name)[i,]
    print(family)
    df = read_delim(paste0(dir,'/mcl/',family,'_mcl_evalue'),  delim = '\t', col_names = FALSE)
    n_clusters = nrow(df)

    # Transpose and reshape dataframe
    df_t = as.data.frame(t(df))
    colnames(df_t) <- paste0("cl_",1:n_clusters) 

    df_longer = df_t %>%
                    melt(variable.names = "cluster", value.name = "protein_id", id.vars=integer()) %>%
                    drop_na() %>%
                    dplyr::rename(cluster = variable)

    # Load FASTA file
    fasta_file <- paste0(dir, '/cd_hit/', family, '.fa') # Replace with the correct FASTA file path
    fasta_sequences <- readAAStringSet(fasta_file) # Load sequences

    # Convert to a named vector
    fasta_data <- as.data.frame(fasta_sequences)
    fasta_df <- tibble(protein_description = names(fasta_sequences), protein_sequence = as.character(fasta_sequences)) %>%
                    mutate(protein_id = str_extract(protein_description, "^[^ ]+"), .before = 1)

    # Match proteins to sequences and add to dataframe
    df_longer <- df_longer %>%
                    left_join(fasta_df, by = "protein_id") %>% tibble()

    # Save the final dataframe with sequences
    write_csv(df_longer, paste0(dir, "/mcl/", family, "_protein_description.csv"))


}

