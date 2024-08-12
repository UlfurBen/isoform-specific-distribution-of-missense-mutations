# Load necessary libraries
library(dplyr)

# Read the BED file
bed_data <- read.table("variants_filtered_sorted.bed", header = FALSE, stringsAsFactors = FALSE, fill=TRUE)

# Read the gene names from the CSV file
gene_names <- read.csv("gene_names.csv", header = FALSE, stringsAsFactors = FALSE)[, 1]

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
  
  # Find all unique values in the entire row that start with "ENST"
  unique_values <- unique(unlist(filtered_data %>% select(everything()) %>% as.vector()))
  unique_values <- unique_values[grepl("^ENST", unique_values)]
  
  # Save the filtered rows in the gene's folder, ensuring tab-separated output
  for (value in unique_values) {
    subset_data <- filtered_data %>% filter_all(any_vars(grepl(value, .)))
    write.table(subset_data, file = file.path(gene_folder, paste0(value, ".bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Print a message for each gene
  cat("Processing complete for gene:", gene, "\n")
}

# Print a final completion message
cat("All processing complete. Filtered data saved in the 'genes' superfolder.\n")

