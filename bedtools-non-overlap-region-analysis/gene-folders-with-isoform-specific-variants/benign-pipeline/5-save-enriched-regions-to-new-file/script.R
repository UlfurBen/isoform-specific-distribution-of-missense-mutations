# Set the paths for the directories and output file
genes_dir <- "~/benign-genes"
results_dir <- "~/results"
output_file <- file.path(results_dir, "benign_filtered_results.bed")

# Create the results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Initialize an empty data frame to collect all filtered rows
all_filtered_data <- data.frame()

# Loop over each folder within the genes directory
folders <- list.dirs(genes_dir, full.names = TRUE, recursive = FALSE)

for (folder in folders) {
  # Define the path to the intersected folder
  intersected_folder <- file.path(folder, "intersected")
  
  # Check if the intersected folder exists
  if (dir.exists(intersected_folder)) {
    # Loop over each file within the intersected folder
    files <- list.files(intersected_folder, full.names = TRUE)
    
    for (file in files) {
      # Read the file
      data <- read.table(file, header = FALSE)
      
      # Filter rows where the value in the 6th column is non-zero
      filtered_data <- data[data[, 6] != 0, ]
      
      # Append the filtered data to the overall data frame
      all_filtered_data <- rbind(all_filtered_data, filtered_data)
    }
  }
}

# Save all filtered rows to a single file in the results directory
write.table(all_filtered_data, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
