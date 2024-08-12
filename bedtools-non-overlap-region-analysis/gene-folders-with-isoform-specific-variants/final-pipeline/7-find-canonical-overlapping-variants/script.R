# Load necessary libraries
library(dplyr)
library(readr)

# Set the superfolder path
superfolder <- "~/genes"

# Load the ENST values for non-canonical isoforms from the CSV file
non_canonical_data <- read_csv("~/60_gene_names_with_enst.csv")

# Get a list of all subfolders within the superfolder
subfolders <- list.dirs(superfolder, full.names = TRUE, recursive = FALSE)

# Function to process isoform files
process_isoforms <- function(subfolder, label, non_canonical_data) {
  # Define the path to the isoform folders
  isoform_folder <- file.path(subfolder, label)
  
  # Get the gene name (subfolder name)
  gene_name <- basename(subfolder)
  
  # Filter the non-canonical ENST values for this gene
  non_canonical_enst <- non_canonical_data %>%
    filter(`X...EM.genes` == gene_name) %>%
    pull(Canonical_ENST)
  
  # Define the paths to the canonical and non-canonical isoform files
  canonical_file <- file.path(isoform_folder, paste0("canonical_", label, ".bed"))
  non_canonical_file <- file.path(isoform_folder, paste0("non_canonical_", label, ".bed"))
  
  # Check if both files exist
  if (file.exists(canonical_file) && file.exists(non_canonical_file)) {
    # Read the isoform files
    canonical_data <- read_tsv(canonical_file, col_names = FALSE)
    non_canonical_data <- read_tsv(non_canonical_file, col_names = FALSE)
    
    # Filter the non-canonical data for the ENST values associated with this gene
    non_canonical_filtered <- non_canonical_data %>%
      filter(X1 %in% non_canonical_enst)  # Assuming ENST values are in the first column
    
    # Find intersections (common rows) with the canonical data
    intersect_data <- semi_join(non_canonical_filtered, canonical_data, by = colnames(non_canonical_filtered))
    
    # Find non-intersecting data
    non_intersect_data <- anti_join(non_canonical_filtered, canonical_data, by = colnames(non_canonical_filtered))
    
    # Create the non-canonical-isoforms folder if it doesn't exist
    output_folder <- file.path(isoform_folder, "non-canonical-isoforms")
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }
    
    # Write the intersected data to a file
    if (nrow(intersect_data) > 0) {
      write_tsv(intersect_data, file.path(output_folder, paste0("intersect_", label, ".bed")), col_names = FALSE)
    }
    
    # Write the non-intersected data to a file
    if (nrow(non_intersect_data) > 0) {
      write_tsv(non_intersect_data, file.path(output_folder, paste0("non_intersect_", label, ".bed")), col_names = FALSE)
    }
    
    # Print a message indicating completion for this label
    cat("Processed", label, "isoforms for subfolder:", gene_name, "\n")
  } else {
    cat("Canonical or non-canonical isoform file not found in subfolder:", gene_name, "for", label, "\n")
  }
}

# Loop over each subfolder
for (subfolder in subfolders) {
  # Process the benign isoforms
  process_isoforms(subfolder, "benign", non_canonical_data)
  
  # Process the pathogenic isoforms
  process_isoforms(subfolder, "pathogenic", non_canonical_data)
}

# Print a final completion message
cat("All processing complete. Intersected and non-intersected files have been created.\n")
