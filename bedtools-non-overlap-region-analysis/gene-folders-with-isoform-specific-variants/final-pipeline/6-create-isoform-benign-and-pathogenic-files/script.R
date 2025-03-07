# Set the superfolder path
superfolder <- "~/genes"

# Get a list of all subfolders within the superfolder
subfolders <- list.dirs(superfolder, full.names = TRUE, recursive = FALSE)

# Function to process benign or pathogenic files
process_files <- function(subfolder, variants_file, label) {
  # Read the variants file
  variants_data <- readLines(variants_file)
  
  # Extract unique ENST values (assumes ENST values are present in the file)
  enst_values <- unique(unlist(regmatches(variants_data, gregexpr("ENST\\d+", variants_data))))
  
  # Create subfolder for the label (benign or pathogenic)
  subfolder_label <- file.path(subfolder, label)
  if (!dir.exists(subfolder_label)) {
    dir.create(subfolder_label)
  }
  
  # Loop over each unique ENST value
  for (enst in enst_values) {
    # Grep for the ENST value in the variants file
    enst_data <- grep(enst, variants_data, value = TRUE)
    
    # Create a file name for the ENST value
    enst_file <- file.path(subfolder_label, paste0(enst, "_", label, ".bed"))
    
    # Write the ENST data to the file
    writeLines(enst_data, enst_file)
    
    # Print a message indicating completion for this ENST value
    cat("Processed ENST:", enst, "in", label, "for subfolder:", basename(subfolder), "\n")
  }
}

# Loop over each subfolder
for (subfolder in subfolders) {
  # Get the folder name
  folder_name <- basename(subfolder)
  
  # Construct the paths to the benign and pathogenic files
  benign_file <- file.path(subfolder, paste0("benign_", folder_name, ".bed"))
  pathogenic_file <- file.path(subfolder, paste0("pathogenic_", folder_name, ".bed"))
  
  # Process the benign file if it exists
  if (file.exists(benign_file)) {
    process_files(subfolder, benign_file, "benign")
  } else {
    cat("No benign file found in subfolder:", folder_name, "\n")
  }
  
  # Process the pathogenic file if it exists
  if (file.exists(pathogenic_file)) {
    process_files(subfolder, pathogenic_file, "pathogenic")
  } else {
    cat("No pathogenic file found in subfolder:", folder_name, "\n")
  }
}

# Print a final completion message
cat("All processing complete. ENST-specific files have been created.\n")
