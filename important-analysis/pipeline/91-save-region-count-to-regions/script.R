# Load necessary library
library(stringr)

# Define the path to the genes folder
genes_folder <- "~/genes"

# Define the output file path for enriched regions in the main ~/genes directory
enriched_regions_file <- file.path(genes_folder, "enriched_regions.bed")

# Initialize an empty data frame to store enriched regions from all genes
all_enriched_regions <- data.frame()

# Get the list of subfolders in the genes folder
subfolders <- list.dirs(genes_folder, recursive = FALSE)

# Loop through each subfolder
for (subfolder in subfolders) {
  
  # Define the path to the region_variant_count folder
  region_variant_count_folder <- file.path(subfolder, "region_variant_count")
  
  # Check if the region_variant_count folder exists
  if (dir.exists(region_variant_count_folder)) {
    
    # Get the list of files in the region_variant_count folder
    count_files <- list.files(region_variant_count_folder, pattern = "\\.txt$", full.names = TRUE)
    
    # Loop through each file in the region_variant_count folder
    for (count_file in count_files) {
      
      # Read the file into a data frame
      variant_data <- read.table(count_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      
      # Append all rows to the all_enriched_regions data frame
      all_enriched_regions <- rbind(all_enriched_regions, variant_data)
    }
  } else {
    cat("Region variant count folder not found in", subfolder, "\n")
  }
}

# Sort the data frame by the 6th column in descending order
all_enriched_regions <- all_enriched_regions[order(-all_enriched_regions$V6), ]

# Write all regions to the enriched_regions.bed file
write.table(all_enriched_regions, enriched_regions_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
cat("Saved all regions, sorted by the 6th column, to", enriched_regions_file, "\n")

# Print a final completion message
cat("Enriched regions have been saved to the main ~/genes directory.\n")
