#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(parallel)

# Set the superfolder path
superfolder <- "~/genes"

# Get a list of all subfolders within the superfolder
subfolders <- list.dirs(superfolder, full.names = TRUE, recursive = FALSE)

# Function to process each folder (benign or pathogenic)
process_folder <- function(folder) {
  # Get a list of BED files in the folder
  bed_files <- list.files(folder, pattern = "\\.bed$", full.names = TRUE)
  
  # Create a new folder for isoform-specific variants
  isoform_specific_folder <- file.path(folder, "isoform_specific_variants")
  if (!dir.exists(isoform_specific_folder)) {
    dir.create(isoform_specific_folder)
  }
  
  # Loop over each BED file
  for (bed_file in bed_files) {
    # Define a temporary file to store non-intersecting variants
    temp_file <- file.path(isoform_specific_folder, paste0("temp_", basename(bed_file)))
    
    # Start with the full content of the file as potentially non-intersecting variants
    system(paste("cp", bed_file, temp_file))
    
    # Loop over all other BED files in the folder
    for (other_file in bed_files) {
      if (other_file != bed_file) {
        # Find non-intersecting variants using bedtools intersect -v
        intersect_command <- paste("bedtools intersect -a", temp_file, "-b", other_file, "-v >", temp_file)
        system(intersect_command)
      }
    }
    
    # Move the final non-intersecting variants to the isoform_specific_variants folder
    final_output <- file.path(isoform_specific_folder, basename(bed_file))
    file.rename(temp_file, final_output)
    
    # Print a message indicating completion for this file
    cat("Processed", basename(bed_file), "and saved isoform-specific variants to", final_output, "\n")
  }
}

# Function to process each subfolder (benign and pathogenic)
process_subfolder <- function(subfolder) {
  # Process the benign folder if it exists
  benign_folder <- file.path(subfolder, "benign")
  if (dir.exists(benign_folder)) {
    cat("Processing benign folder in", subfolder, "\n")
    process_folder(benign_folder)
  }
  
  # Process the pathogenic folder if it exists
  pathogenic_folder <- file.path(subfolder, "pathogenic")
  if (dir.exists(pathogenic_folder)) {
    cat("Processing pathogenic folder in", subfolder, "\n")
    process_folder(pathogenic_folder)
  }
}

# Set up the cluster for parallel processing
num_cores <- detectCores() - 1  # Use all available cores except one
cl <- makeCluster(num_cores)

# Load necessary libraries on each worker node
clusterEvalQ(cl, {
  library(data.table)
})

# Export necessary variables and functions to the cluster
clusterExport(cl, varlist = c("process_folder", "process_subfolder", "system"))

# Apply the function in parallel to each subfolder
parLapply(cl, subfolders, process_subfolder)

# Stop the cluster
stopCluster(cl)

# Print a final completion message
cat("All non-intersecting variant processing is complete. Isoform-specific variants have been saved.\n")
