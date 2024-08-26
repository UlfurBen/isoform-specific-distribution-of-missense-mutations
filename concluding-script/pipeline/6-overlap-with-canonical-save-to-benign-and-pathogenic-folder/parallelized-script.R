#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(parallel)

# Set the superfolder path
superfolder <- "~/genes"

# Load the important genes with canonical ENST values
important_genes <- fread("~/important_genes_with_enst.csv")

# Get a list of all subfolders within the superfolder
subfolders <- list.dirs(superfolder, full.names = TRUE, recursive = FALSE)

# Function to process each subfolder
process_subfolder <- function(subfolder) {
  # Extract the gene name from the subfolder name
  gene_name <- basename(subfolder)
  
  # Check if the gene name is in the important genes list
  gene_info <- important_genes[important_genes[[1]] == gene_name]
  
  if (nrow(gene_info) == 1) {
    # Get the canonical ENST value for this gene
    canonical_enst <- gene_info[[2]]
    
    # Define the path to the canonical.bed file
    canonical_file <- file.path(subfolder, "canonical.bed")
    
    # Check if the canonical.bed file exists
    if (file.exists(canonical_file)) {
      
      # Create a new folder within the subfolder for the intersection results
      output_folder <- file.path(subfolder, "intersections")
      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
      }
      
      # Create subfolders for benign and pathogenic results
      benign_folder <- file.path(subfolder, "benign")
      pathogenic_folder <- file.path(subfolder, "pathogenic")
      if (!dir.exists(benign_folder)) {
        dir.create(benign_folder)
      }
      if (!dir.exists(pathogenic_folder)) {
        dir.create(pathogenic_folder)
      }
      
      # Get a list of non-canonical BED files in the subfolder
      bed_files <- list.files(subfolder, pattern = "\\.bed$", full.names = TRUE)
      
      # Filter out the canonical.bed and variants_Crebbp.bed files from the list
      bed_files <- bed_files[!(basename(bed_files) %in% c("canonical.bed", paste0("variants_", gene_name, ".bed")))]
      
      # Loop over each non-canonical BED file
      for (bed_file in bed_files) {
        # Define the output file name for intersection
        intersect_file <- file.path(output_folder, paste0(basename(bed_file), "_labeled.bed"))
        
        # Construct the bedtools intersect command
        intersect_command <- paste("bedtools intersect -a", canonical_file, "-b", bed_file, "-wa >", intersect_file)
        
        # Run the bedtools intersect command
        system(intersect_command)
        
        # Print a message indicating completion for this intersection
        cat("Intersected", basename(bed_file), "with canonical.bed. Output saved to", intersect_file, "\n")
        
        # Grep case-insensitive 'benign' excluding lines with 'pathogenic' in the intersected data
        benign_output <- file.path(benign_folder, paste0(basename(bed_file), "_benign.bed"))
        grep_benign_command <- paste("grep -i 'benign' ", intersect_file, " | grep -vi 'pathogenic' >", benign_output)
        system(grep_benign_command)
        cat("Grep 'benign' (excluding 'pathogenic') in", intersect_file, ". Output saved to", benign_output, "\n")
        
        # Grep case-insensitive 'pathogenic' excluding lines with 'benign' in the intersected data
        pathogenic_output <- file.path(pathogenic_folder, paste0(basename(bed_file), "_pathogenic.bed"))
        grep_pathogenic_command <- paste("grep -i 'pathogenic' ", intersect_file, " | grep -vi 'benign' >", pathogenic_output)
        system(grep_pathogenic_command)
        cat("Grep 'pathogenic' (excluding 'benign') in", intersect_file, ". Output saved to", pathogenic_output, "\n")
      }
    } else {
      cat("Canonical file not found for gene:", gene_name, "\n")
    }
  } else {
    cat("No canonical ENST found for gene:", gene_name, "\n")
  }
}

# Set up the cluster for parallel processing
num_cores <- 30
cl <- makeCluster(num_cores)

# Load necessary libraries on each worker node
clusterEvalQ(cl, {
  library(data.table)
})

# Export necessary variables and functions to the cluster
clusterExport(cl, varlist = c("important_genes", "process_subfolder", "superfolder", "fread", "system"))

# Apply the function in parallel to each subfolder
parLapply(cl, subfolders, process_subfolder)

# Stop the cluster
stopCluster(cl)

# Print a final completion message
cat("All intersections and grep operations are complete. Output files have been created in respective subfolders.\n")
