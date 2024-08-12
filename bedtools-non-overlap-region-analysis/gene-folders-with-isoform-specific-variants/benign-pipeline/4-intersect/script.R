# Load necessary libraries
library(dplyr)

# Define the path to the genes directory
genes_dir <- "~/benign-genes"

# Function to perform intersection
intersect_files <- function(regions_file, uniq_file, output_file) {
  if (file.info(uniq_file)$size > 0) {  # Check if uniq_file is not empty
    # Perform intersection using bedtools
    system(paste("bedtools intersect -c -a", regions_file, "-b", uniq_file, ">", output_file))
  } else {
    cat("Skipping intersection for", basename(regions_file), "- unique_intervals file is empty.\n")
  }
}

# Loop over each subfolder in genes_dir
for (subfolder in list.dirs(genes_dir, full.names = TRUE, recursive = FALSE)) {
  regions_dir <- file.path(subfolder, "regions")
  unique_intervals_dir <- file.path(subfolder, "unique_intervals")
  
  # Create intersected folder if it doesn't exist
  intersected_dir <- file.path(subfolder, "intersected")
  if (!dir.exists(intersected_dir)) {
    dir.create(intersected_dir)
  }
  
  # Get list of region files
  region_files <- list.files(regions_dir, pattern = "_regions.bed$", full.names = TRUE)
  
  # Process each regions file
  for (region_file in region_files) {
    # Extract the base name to match with unique_intervals file
    base_name <- gsub("_regions.bed$", "", basename(region_file))
    
    # Construct the unique_intervals file path
    uniq_file <- file.path(unique_intervals_dir, paste0(base_name, "_uniq.bed"))
    
    # Check if the corresponding unique_intervals file exists
    if (file.exists(uniq_file)) {
      # Construct the output file path
      output_file <- file.path(intersected_dir, paste0(base_name, "_intersect.bed"))
      
      # Perform the intersection
      intersect_files(region_file, uniq_file, output_file)
    } else {
      cat("Unique intervals file does not exist for", base_name, "\n")
    }
  }
}

cat("Intersection process completed.\n")
