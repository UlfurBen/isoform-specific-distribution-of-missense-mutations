# Load necessary libraries
library(dplyr)

# Define the path to the superfolder containing the genes
superfolder <- "~/genes"

# Define the path to the variants.bed file
variants_file <- file.path(superfolder, "variants.bed")

# Read the names of the folders within the superfolder (assumes folder names are the gene names)
gene_folders <- list.dirs(superfolder, full.names = FALSE, recursive = FALSE)

# Read the variants.bed file
variants_data <- read.table(variants_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# Loop over each gene folder name
for (gene in gene_folders) {
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

# Print a final completion message
cat("All processing complete. Filtered data saved in the respective gene folders.\n")
