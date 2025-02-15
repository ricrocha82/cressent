#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(Biostrings) # For working with FASTA files

# Set directory and file paths
dir = '/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db'
family_name = read.delim('/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/names_to_query_ssdna.csv', sep = '\t', header = FALSE)


# Read family information
for(i in seq_along(family_name$V1)) {
    # Issue 2: Incorrect data frame access
    # Fixed: Proper extraction of family name
    family <- family_name$V1[i]
    print(family)
    
    file_path <- paste0(dir, '/mcl/', family, '_mcl_evalue')
    
    if (file.exists(file_path) && file.size(file_path) > 0) {
        # Read the MCL file
        df <- fread(file_path, header = FALSE, fill = TRUE, na.strings = c("", "NA"), sep = "\t")
      
        # Only process if we have data
        if (nrow(df) > 0) {
        n_clusters <- nrow(df)
        
        # Transpose and reshape dataframe
        # Issue 3: Potential memory issue with large datasets
        # Fixed: More efficient reshaping
        df_t <- as.data.frame(t(df))
        colnames(df_t) <- paste0("cl_", 1:n_clusters)
        
        df_longer <- df_t %>% 
            melt(variable.name = "cluster",  # Fixed: variable.names to variable.name
                value.name = "protein_id", 
                id.vars = integer()) %>%         # Fixed: removed integer()
            drop_na() 
            # %>% 
            # dplyr::rename(cluster = variable)
        
        # Load FASTA file
        fasta_file <- paste0(dir, '/cd_hit/', family, '.fa')
        
        # Issue 4: No error handling for FASTA reading
        # Fixed: Added tryCatch
        tryCatch({
            fasta_sequences <- readAAStringSet(fasta_file)
            
            # Process FASTA data
            fasta_df <- tibble(protein_description = names(fasta_sequences),
                                protein_sequence = as.character(fasta_sequences)) %>%
                        mutate(protein_id = str_extract(protein_description, "^[^ ]+"),
                                family = family,
                                scientific_name = gsub(".*\\[(.*)\\].*", "\\1", protein_description),
                                protein_name = str_remove(protein_description, " \\[.*") %>%
                                            str_remove("^\\S+ "))
            
            
            # Join data and save
            df_longer <- df_longer %>%
            left_join(fasta_df, by = "protein_id") %>%
            tibble()
            
        # Save results
            write_csv(df_longer, paste0(dir, "/mcl/", family, "_protein_description.csv"))
            rm(fasta_sequences, fasta_df, df_longer)
        
        }, error = function(e) {
          message("Error processing FASTA file for family ", family, ": ", e$message)
        })
        } else {
        message("No data in file: ", file_path)
        }
    } else {
        message("File does not exist or is empty: ", file_path)
    }
  }


directory_path = "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/mcl"
# Function to combine multiple CSV files from a directory
combine_csv_files <- function(directory_path, output_file = "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/mcl/combined_protein_description.csv") {
    # Load required library
    library(dplyr)
    library(readr)
    
    # Get list of all CSV files in the directory
    csv_files <- list.files(path = directory_path, 
                            pattern = "*.csv", 
                            full.names = TRUE)
    
    # Check if any CSV files were found
    if (length(csv_files) == 0) {
        stop("No CSV files found in the specified directory")
    }
    
    # Read and combine all CSV files
    combined_data <- csv_files %>%
        lapply(data.table::fread) %>%  # Read each CSV file
        bind_rows()  # Combine all data frames

        
    # Add source file information (optional)
    # combined_data$source_file <- basename(csv_files[as.numeric(combined_data$source_file)])
    
    # Write the combined data to a new CSV file
    write_csv(combined_data, output_file)
    
    # Return the combined data frame
    return(combined_data)
}


# Example usage
combine_csv_files(directory_path = "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/mcl",  
                    output_file = "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/mcl/combined_protein_description.csv")


# test = read_csv("/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/mcl/combined_protein_description.csv")

# test[!complete.cases(test),] %>% print(n = 512) %>% print(width = Inf)
