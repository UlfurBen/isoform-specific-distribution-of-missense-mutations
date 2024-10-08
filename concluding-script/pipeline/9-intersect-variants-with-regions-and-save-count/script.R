# Load necessary library
library(stringr)

# Define the path to the genes folder
genes_folder <- "~/genes"

# Get the list of subfolders in the genes folder
subfolders <- list.dirs(genes_folder, recursive = FALSE)

# Loop through each subfolder
for (subfolder in subfolders) {
  
  # Get the gene name from the subfolder name
  gene_name <- basename(subfolder)
  
  # Define paths to the regions, benign, and pathogenic folders
  regions_folder <- file.path(subfolder, "regions")
  benign_folder <- file.path(subfolder, "benign", "isoform_specific_variants")
  pathogenic_folder <- file.path(subfolder, "pathogenic", "isoform_specific_variants")
  
  # Create a folder for region_variant_count if it doesn't exist
  region_variant_count_folder <- file.path(subfolder, "region_variant_count")
  if (!dir.exists(region_variant_count_folder)) {
    dir.create(region_variant_count_folder)
  }
  
  # Check if the regions folder exists
  if (dir.exists(regions_folder)) {
    
    # Get the list of .bed files in the regions folder
    region_files <- list.files(regions_folder, pattern = "\\_regions\\.bed$", full.names = TRUE)
    
    # Loop through each .bed file in the regions folder
    for (region_file in region_files) {
      
      # Extract the ENST value from the file name
      enst_value <- str_extract(basename(region_file), "ENST\\d+")
      
      # Define the corresponding pathogenic and benign bed file paths
      pathogenic_file <- file.path(pathogenic_folder, paste0(enst_value, "_pathogenic.bed"))
      benign_file <- file.path(benign_folder, paste0(enst_value, "_benign.bed"))
      
      # Define the output file path for the combined variant count results
      combined_count_output <- file.path(region_variant_count_folder, paste0(enst_value, "_region_variant_count.txt"))
      
      # Initialize a data frame to store region coordinates and variant counts
      region_data <- read.table(region_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      region_data$Benign_Count <- 0
      region_data$Pathogenic_Count <- 0
      
      # Count the number of benign variants in each region if the file exists
      if (file.exists(benign_file)) {
        benign_counts <- system(paste("bedtools intersect -a", region_file, "-b", benign_file, "-c"), intern = TRUE)
        benign_counts <- read.table(text = benign_counts, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        region_data$Benign_Count <- benign_counts$V4
        cat("Calculated benign variant counts for", region_file, "\n")
      } else {
        cat("Benign file not found for ENST:", enst_value, "in", subfolder, "\n")
      }
      
      # Count the number of pathogenic variants in each region if the file exists
      if (file.exists(pathogenic_file)) {
        pathogenic_counts <- system(paste("bedtools intersect -a", region_file, "-b", pathogenic_file, "-c"), intern = TRUE)
        pathogenic_counts <- read.table(text = pathogenic_counts, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        region_data$Pathogenic_Count <- pathogenic_counts$V4
        cat("Calculated pathogenic variant counts for", region_file, "\n")
      } else {
        cat("Pathogenic file not found for ENST:", enst_value, "in", subfolder, "\n")
      }
      
      # Write the combined counts to the output file
      write.table(region_data, combined_count_output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      cat("Saved combined variant counts to", combined_count_output, "\n")
    }
  } else {
    cat("Regions folder not found in", subfolder, "\n")
  }
}

# Print a final completion message
cat("All variant counting is complete. Results have been saved in the respective region_variant_count folders.\n")
