# Load necessary libraries
library(stringr)
library(parallel)

# Define the path to the genes folder and the Homo_sapiens.GRCh37.87.bed file
genes_folder <- "~/genes"
homo_sapiens_file <- "filtered_no_scientific_notation.bed"

# Get the list of subfolders in the genes folder
subfolders <- list.dirs(genes_folder, recursive = FALSE)

# Function to process each subfolder
process_subfolder <- function(subfolder) {
  # Get the gene name from the subfolder name
  gene_name <- basename(subfolder)
  
  # Get the list of .bed files in the subfolder
  bed_files <- list.files(subfolder, pattern = "\\.bed$", full.names = TRUE)
  
  # Loop through each .bed file in the subfolder
  for (bed_file in bed_files) {
    
    # Get the base name of the bed file without extension
    bed_file_name <- tools::file_path_sans_ext(basename(bed_file))
    
    # Skip the file if it is canonical.bed or variants_{gene_name}.bed
    if (bed_file_name %in% c("canonical", paste0("variants_", gene_name))) {
      next
    }
    
    # Grep the gene name and the bed file name from the Homo_sapiens.GRCh37.87.bed file
    grep_command <- paste0("grep -i ", gene_name, " ", homo_sapiens_file, " | grep ", bed_file_name)
    grep_results <- system(grep_command, intern = TRUE)
    
    # If there are results, write them to a new file in the regions folder
    if (length(grep_results) > 0) {
      # Create the regions folder if it doesn't exist
      regions_folder <- file.path(subfolder, "regions")
      if (!dir.exists(regions_folder)) {
        dir.create(regions_folder)
      }
      
      # Define the output file path
      output_file <- file.path(regions_folder, paste0(bed_file_name, "_regions.bed"))
      
      # Write the grep results to the output file
      writeLines(grep_results, output_file)
    }
  }
}

# Set up the cluster for parallel processing
num_cores <- 30
cl <- makeCluster(num_cores)

# Load necessary libraries on each worker node
clusterEvalQ(cl, {
  library(stringr)
})

# Export necessary variables and functions to the cluster
clusterExport(cl, varlist = c("genes_folder", "homo_sapiens_file", "process_subfolder", "system"))

# Apply the function in parallel to each subfolder
parLapply(cl, subfolders, process_subfolder)

# Stop the cluster
stopCluster(cl)

# Print a final completion message
cat("All grep operations are complete. Regions files have been saved in respective subfolders.\n")
