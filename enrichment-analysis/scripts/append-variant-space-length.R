# Read the enrichment.txt file
enrichment <- read.csv("enrichment.txt", stringsAsFactors = FALSE)

# Add the new column 'variant_space' to the enrichment dataframe
enrichment$variant_space <- NA

# Function to calculate the difference between the last and first line values
calculate_difference <- function(file_path) {
  # Read the file
  data <- read.table(file_path, header = FALSE)
  
  # Check if the file has at least two lines
  if (nrow(data) < 2) {
    cat("File", file_path, "does not have enough data.\n")
    return(NA)
  }
  
  # Calculate the difference
  difference <- data[nrow(data), 1] - data[1, 1]
  
  return(difference)
}

# Directory containing the identifier files
directory <- "1-identifier-genomic-location"

# Loop through each identifier in the enrichment file and process it
for (i in 1:nrow(enrichment)) {
  identifier <- enrichment$Identifier[i]
  file_path <- file.path(directory, paste0(identifier, ".txt"))
  
  # Check if the file exists
  if (file.exists(file_path)) {
    cat("Processing identifier:", identifier, "\n")
    difference <- calculate_difference(file_path)
    
    # Check if the difference is not NA
    if (!is.na(difference)) {
      # Append the difference to the 'variant_space' column
      enrichment$variant_space[i] <- difference
    }
  } else {
    cat("File for identifier", identifier, "not found.\n")
  }
}

# Save the updated enrichment data to enrichment.txt
write.table(enrichment, "enrichment.txt", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Updated enrichment.txt successfully.\n")
