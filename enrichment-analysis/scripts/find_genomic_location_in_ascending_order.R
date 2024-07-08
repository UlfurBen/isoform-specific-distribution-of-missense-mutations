# Load necessary library
if(!dir.exists("1-identifier-genomic-location")){
  dir.create("1-identifier-genomic-location")
}

# Read the enrichment.txt file
enrichment <- read.csv("enrichment.txt", stringsAsFactors = FALSE)

# Function to process the homo_sapiens_variation_missense_EM_genes_phenotype.txt for a given identifier
process_identifier <- function(identifier, em_genes_data) {
  # Filter rows based on the identifier in the second column
  filtered_data <- em_genes_data[grep(identifier, em_genes_data[, 2]), ]
  
  # Check if filtered_data is empty
  if (nrow(filtered_data) == 0) {
    cat("No data found for identifier:", identifier, "\n")
    return(NULL)
  }
  
  cat("Found data for identifier:", identifier, "\n")
  
  # Ensure the 10th column exists and has valid data
  if (ncol(filtered_data) < 10) {
    cat("The data for identifier:", identifier, "does not have 10 columns.\n")
    return(NULL)
  }
  
  # Extract the 10th column
  positions <- filtered_data[, 10]
  
  # Check if positions column is empty or has invalid values
  valid_positions <- grep("^NC_.*g\\..*>", positions)
  if (length(valid_positions) == 0) {
    cat("No valid positions found for identifier:", identifier, "\n")
    return(NULL)
  }
  
  # Filter only valid positions
  positions <- positions[valid_positions]
  
  # Print the positions for debugging
  cat("Positions for identifier", identifier, ":", positions, "\n")
  
  # Strip the value to only contain the numeric part of the genomic location
  stripped_values <- as.numeric(sub(".*g\\.([0-9]+)[A-Z]>[A-Z].*", "\\1", positions))
  
  # Print the stripped values for debugging
  cat("Stripped values for identifier", identifier, ":", stripped_values, "\n")
  
  # Check if stripped_values is empty or contains NA
  if (length(stripped_values) == 0 || all(is.na(stripped_values))) {
    cat("No valid stripped values for identifier:", identifier, "\n")
    return(NULL)
  }
  
  # Sort the values in ascending order
  sorted_values <- sort(stripped_values, na.last = NA)
  
  return(sorted_values)
}

# Read the homo_sapiens_variation_missense_EM_genes_phenotype.txt file
em_genes_data <- read.table("homo_sapiens_variation_missense_EM_genes_phenotype.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# Print first few rows for debugging
cat("First few rows of the EM genes phenotype data:\n")
print(head(em_genes_data))

# Loop through each identifier in the enrichment file and process it
for (identifier in enrichment$Identifier) {
  cat("Processing identifier:", identifier, "\n")
  
  sorted_values <- process_identifier(identifier, em_genes_data)
  
  # Check if sorted_values is not NULL
  if (!is.null(sorted_values)) {
    # Save the result to a file named after the identifier in the new directory
    output_file <- paste0("1-identifier-genomic-location/", identifier, ".txt")
    write.table(sorted_values, file = output_file, row.names = FALSE, col.names = FALSE)
    cat("Saved sorted values for identifier:", identifier, "to file:", output_file, "\n")
  } else {
    cat("No data to save for identifier:", identifier, "\n")
  }
}
