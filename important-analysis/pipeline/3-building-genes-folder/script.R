# Load necessary libraries
library(dplyr)

# Read the BED file
bed_data <- read.table("variants_benign_pathogenic_non_vus_non_conflicting_including_non_annotation.bed", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# Read the gene names from the CSV file (assuming it contains a header)
gene_names <- read.csv("important_genes_with_enst.csv", header = TRUE, stringsAsFactors = FALSE)[, 1]

# Create the superfolder "genes"
superfolder <- "genes"
if (!dir.exists(superfolder)) {
  dir.create(superfolder)
}

# Initialize an empty data frame to store the filtered variants
all_filtered_data <- data.frame()

# Loop over each gene name
for (gene in gene_names) {
  # Filter the data for the current gene name (case insensitive) in the fifth column
  filtered_data <- bed_data %>%
    filter(grepl(gene, V5, ignore.case = TRUE))
  
  # Append the filtered data to the combined data frame
  all_filtered_data <- rbind(all_filtered_data, filtered_data)
  
  # Create an empty folder for the current gene inside the superfolder
  gene_folder <- file.path(superfolder, gene)
  if (!dir.exists(gene_folder)) {
    dir.create(gene_folder)
  }
  
  # Print a message for each gene
  cat("Folder created for gene:", gene, "\n")
}

# Define the output file path for the combined filtered data
output_file <- file.path(superfolder, "variants.bed")

# Write the combined filtered data to a file
write.table(all_filtered_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Print a final completion message
cat("All processing complete. Filtered data saved in 'variants.bed' within the 'genes' superfolder.\n")
