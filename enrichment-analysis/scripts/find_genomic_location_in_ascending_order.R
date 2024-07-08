# Read the enrichment.txt file
enrichment <- read.csv("enrichment.txt")

# Function to process the homo_sapiens_variation_missense_EM_genes.txt for a given identifier
process_identifier <- function(identifier, em_genes_data) {
  # Filter rows based on the identifier in the second column
  filtered_data <- em_genes_data[grep(identifier, em_genes_data[, 2]), ]
  
  # Check if filtered_data is empty
  if (nrow(filtered_data) == 0) {
    cat("No data found for identifier:", identifier, "\n")
    return(NULL)
  }
  
  # Extract the 10th column
  positions <- filtered_data[, 10]
  
  # Check if positions column is empty
  if (length(positions) == 0) {
    cat("No positions found for identifier:", identifier, "\n")
    return(NULL)
  }
  
  # Strip the value to only contain the numeric part
  stripped_values <- as.numeric(sub(".*g\\.([0-9]+)[A-Z]>[A-Z].*", "\\1", positions))
  
  # Check if stripped_values is empty or contains NA
  if (length(stripped_values) == 0 || all(is.na(stripped_values))) {
    cat("No valid stripped values for identifier:", identifier, "\n")
    return(NULL)
  }
  
  # Sort the values in ascending order
  sorted_values <- sort(stripped_values, na.last = NA)
  
  return(sorted_values)
}

# Read the homo_sapiens_variation_missense_EM_genes.txt file
em_genes_data <- read.table("homo_sapiens_variation_missense_EM_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Loop through each identifier in the enrichment file and process it
for (identifier in enrichment$Identifier) {
  cat("Processing identifier:", identifier, "\n")
  
  sorted_values <- process_identifier(identifier, em_genes_data)
  
  # Check if sorted_values is not NULL
  if (!is.null(sorted_values)) {
    # Save the result to a file named after the identifier
    write.table(sorted_values, file = paste0(identifier, ".txt"), row.names = FALSE, col.names = FALSE)
    cat("Saved sorted values for identifier:", identifier, "\n")
  } else {
    cat("No data to save for identifier:", identifier, "\n")
  }
}
