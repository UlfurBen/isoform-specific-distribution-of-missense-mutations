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

# Function to process variants files
process_files <- function(variants_file, canonical_enst, canonical_data_file, subfolder) {
  # Read the variants file
  variants_data <- readLines(variants_file)
  
  # Extract unique ENST values (assumes ENST values are present in the file)
  enst_values <- unique(unlist(regmatches(variants_data, gregexpr("ENST\\d+", variants_data))))
  
  # Loop over each unique ENST value
  for (enst in enst_values) {
    if (enst == canonical_enst) {
      # If ENST is canonical, read from canonical data file
      canonical_variants <- readLines(canonical_data_file)
      enst_data <- grep(enst, canonical_variants, value = TRUE)
      
      # Save to a file named canonical.bed
      enst_file <- file.path(subfolder, "canonical.bed")
    } else {
      # Otherwise, grep from the current variants file
      enst_data <- grep(enst, variants_data, value = TRUE)
      
      # Create a file name for the ENST value
      enst_file <- file.path(dirname(variants_file), paste0(enst, ".bed"))
    }
    
    # Write the ENST data to the file
    writeLines(enst_data, enst_file)
    
    # Print a message indicating completion for this ENST value
    cat("Processed ENST:", enst, "in file:", basename(enst_file), "\n")
  }
}

# Function to process each subfolder in parallel
process_subfolder <- function(subfolder) {
  # Extract the gene name from the subfolder name
  gene_name <- basename(subfolder)
  
  # Check if the gene name is in the important genes list
  gene_info <- important_genes[important_genes[[1]] == gene_name]
  
  if (nrow(gene_info) == 1) {
    # Get the canonical ENST value for this gene
    canonical_enst <- gene_info[[2]]
    
    # Get the list of BED files in the subfolder
    bed_files <- list.files(subfolder, pattern = "\\.bed$", full.names = TRUE)
    
    # Process each BED file
    for (bed_file in bed_files) {
      # Use the canonical data file for the canonical ENST
      canonical_data_file <- "~/variants_benign_pathogenic_non_vus_non_conflicting.bed"
      process_files(bed_file, canonical_enst, canonical_data_file, subfolder)
    }
  } else {
    cat("No canonical ENST found for gene:", gene_name, "\n")
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
clusterExport(cl, varlist = c("important_genes", "process_files", "process_subfolder", "superfolder"))

# Apply the function in parallel to each subfolder
parLapply(cl, subfolders, process_subfolder)

# Stop the cluster
stopCluster(cl)

# Print a final completion message
cat("All processing complete. ENST-specific files have been created.\n")
