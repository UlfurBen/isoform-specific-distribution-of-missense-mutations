# Load necessary libraries
library(dplyr)
library(parallel)

# Define the path to the superfolder containing the genes
superfolder <- "~/genes"

# Define the path to the variants.bed file
variants_file <- file.path(superfolder, "variants.bed")

# Read the names of the folders within the superfolder (assumes folder names are the gene names)
gene_folders <- list.dirs(superfolder, full.names = FALSE, recursive = FALSE)

# Read the variants.bed file
variants_data <- read.table(variants_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# Function to process each gene folder
process_gene <- function(gene) {
  # Filter the data for the current gene name (case insensitive) in the fifth column
  filtered_data <- variants_data %>%
    filter(grepl(gene, V5, ignore.case = TRUE))
  
  # Define the output file path for the filtered data
  output_file <- file.path(superfolder, gene, paste0("variants_", gene, ".bed"))
  
  # Write the filtered data to a file
  write.table(filtered_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Print a message for each gene
  cat("Variants for gene:", gene, "saved to", output_file, "\n")
}

# Set up the cluster for parallel processing
num_cores <- detectCores() - 1  # Use all available cores except one
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, varlist = c("variants_data", "superfolder", "process_gene"))

# Apply the function in parallel to each gene folder
parLapply(cl, gene_folders, process_gene)

# Stop the cluster
stopCluster(cl)

# Print a final completion message
cat("All processing complete. Filtered data saved in the respective gene folders.\n")
