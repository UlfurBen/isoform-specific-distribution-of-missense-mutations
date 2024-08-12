# Load necessary libraries
library(dplyr)

# Read the BED file
bed_data <- read.table("variants_benign_pathogenic_non_vus_non_conflicting.bed", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# Read the gene names from the CSV file
gene_names <- read.csv("60_EM_gene_names.csv", header = FALSE, stringsAsFactors = FALSE)[, 1]

# Create the superfolder "genes"
superfolder <- "genes"
if (!dir.exists(superfolder)) {
  dir.create(superfolder)
}

# Loop over each gene name
for (gene in gene_names) {
  # Filter the data for the current gene name (case insensitive) in the fifth column
  filtered_data <- bed_data %>%
    filter(grepl(gene, V5, ignore.case = TRUE))
  
  # Create a folder for the current gene inside the superfolder
  gene_folder <- file.path(superfolder, gene)
  if (!dir.exists(gene_folder)) {
    dir.create(gene_folder)
  }
  
  # Define the output file path for the filtered data
  output_file <- file.path(gene_folder, paste0(gene, "_filtered.bed"))
  
  # Write the filtered data to a file
  write.table(filtered_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Print a message for each gene
  cat("Processing complete for gene:", gene, "\n")
}

# Print a final completion message
cat("All processing complete. Filtered data saved in the 'genes' superfolder.\n")
