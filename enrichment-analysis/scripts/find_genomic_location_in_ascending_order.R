# Read the enrichment.txt file
enrichment <- read.csv("enrichment.txt")

# Function to process the homo_sapiens_variation_missense_ClinVar.txt for a given identifier
process_identifier <- function(identifier, clinvar_data) {
  # Filter rows based on the identifier
  filtered_data <- clinvar_data[grep(identifier, clinvar_data[, 1]), ]
  
  # Extract the 10th column
  positions <- filtered_data[, 10]
  
  # Strip the value to only contain the numeric part
  stripped_values <- as.numeric(sub(".*g\\.([0-9]+)[A-Z]>[A-Z].*", "\\1", positions))
  
  # Sort the values in ascending order
  sorted_values <- sort(stripped_values)
  
  return(sorted_values)
}

# Read the homo_sapiens_variation_missense_ClinVar.txt file
clinvar_data <- read.table("homo_sapiens_variation_missense_ClinVar.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Loop through each identifier in the enrichment file and process it
for (identifier in enrichment$Identifier) {
  sorted_values <- process_identifier(identifier, clinvar_data)
  
  # Save the result to a file named after the identifier
  write.table(sorted_values, file = paste0(identifier, ".txt"), row.names = FALSE, col.names = FALSE)
}
