# Load necessary library
library(dplyr)

# Define the function to sort and intersect BED files
process_bed_files <- function(input_folder, output_folder) {
  # Create the output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }

  # List all subfolders in the input folder
  subfolders <- list.dirs(input_folder, recursive = FALSE)
  
  # Process each subfolder
  for (subfolder in subfolders) {
    # Get the ENST value (subfolder name)
    enst_value <- basename(subfolder)
    
    # Create corresponding subfolder in the output folder
    output_subfolder <- file.path(output_folder, enst_value)
    if (!dir.exists(output_subfolder)) {
      dir.create(output_subfolder)
    }
    
    # List all BED files in the current subfolder
    bed_files <- list.files(subfolder, pattern = "\\.bed$", full.names = TRUE)
    
    # Sort each BED file using bedtools
    for (bed_file in bed_files) {
      sorted_bed_file <- bed_file
      system(paste("bedtools sort -i", bed_file, ">", sorted_bed_file))
    }
    
    # Find unique variants in each BED file
    for (bed_file in bed_files) {
      other_bed_files <- setdiff(bed_files, bed_file)
      unique_variants_file <- file.path(output_subfolder, basename(bed_file))
      
      # Intersect with other BED files to find unique variants
      intersect_command <- paste("bedtools intersect -v -a", bed_file, "-b", paste(other_bed_files, collapse = " "), ">", unique_variants_file)
      system(intersect_command)
    }
  }
}

# Define the input and output folders
input_folder <- "crebbp"
output_folder <- "crebbp-isoform-specific"

# Process the BED files
process_bed_files(input_folder, output_folder)

# Output message to confirm completion
cat("Processing complete. Unique variants saved in", output_folder, "\n")
