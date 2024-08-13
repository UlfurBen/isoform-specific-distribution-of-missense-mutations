#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)

# Set the superfolder path
superfolder <- "~/genes"

# Get a list of all subfolders within the superfolder
subfolders <- list.dirs(superfolder, full.names = TRUE, recursive = FALSE)

# Function to process variants files
process_files <- function(variants_file) {
  # Read the variants file
  variants_data <- readLines(variants_file)
  
  # Extract unique ENST values (assumes ENST values are present in the file)
  enst_values <- unique(unlist(regmatches(variants_data, gregexpr("ENST\\d+", variants_data))))
  
  # Loop over each unique ENST value
  for (enst in enst_values) {
    # Grep for the ENST value in the variants file
    enst_data <- grep(enst, variants_data, value = TRUE)
    
    # Create a file name for the ENST value
    enst_file <- file.path(dirname(variants_file), paste0(enst, ".bed"))
    
    # Write the ENST data to the file
    writeLines(enst_data, enst_file)
    
    # Print a message indicating completion for this ENST value
    cat("Processed ENST:", enst, "in file:", basename(variants_file), "\n")
  }
}

# Loop over each subfolder
for (subfolder in subfolders) {
  # Get a list of BED files in the subfolder
  bed_files <- list.files(subfolder, pattern = "\\.bed$", full.names = TRUE)
  
  # Process each BED file
  for (bed_file in bed_files) {
    process_files(bed_file)
  }
}

# Print a final completion message
cat("All processing complete. ENST-specific files have been created.\n")
