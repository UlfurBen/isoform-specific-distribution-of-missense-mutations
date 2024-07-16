# Define file paths
gene_csv_file <- "The-Epigenetic-Machinery.csv"
variation_file <- "homo_sapiens_variation_missense_ClinVar.txt"
output_file <- "count-ranges.txt"

# Load gene names from a CSV file
gene_data <- read.csv(gene_csv_file, header = TRUE)

# Initialize a data frame to store the counts for each gene
range_counts_df <- data.frame(Gene = character(), IdentifierType = character(), Range = character(), Count = integer(), stringsAsFactors = FALSE)

# Define the count ranges
ranges <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-150", "151-200", "201-1000", ">1000")

# Process each gene name from the CSV file
for (gene_name in gene_data[, 1]) {
  
  # Extract identifiers related to the current gene from the variation file
  command <- paste0("grep -i -w '", gene_name, "' ", variation_file, " | awk -F '\t' '{print $23}' | sort -u > temp_identifiers.txt")
  system(command)
  
  # Read the identifiers from the temporary file
  temp_identifiers <- readLines("temp_identifiers.txt")
  
  # Separate canonical and non-canonical identifiers
  canonical_identifiers <- temp_identifiers[!grepl("-", temp_identifiers)]
  non_canonical_identifiers <- temp_identifiers[grepl("-", temp_identifiers)]
  
  # Function to count occurrences of identifiers in the variation file
  count_occurrences <- function(identifiers, type) {
    range_counts <- setNames(rep(0, length(ranges)), ranges)
    
    for (id in identifiers) {
      command <- paste0("grep -c -w '", id, "' ", variation_file)
      count <- as.integer(system(command, intern = TRUE))
      
      # Categorize the counts into the predefined ranges
      if (count >= 1 && count <= 5) {
        range_counts["1-5"] <- range_counts["1-5"] + 1
      } else if (count >= 6 && count <= 10) {
        range_counts["6-10"] <- range_counts["6-10"] + 1
      } else if (count >= 11 && count <= 20) {
        range_counts["11-20"] <- range_counts["11-20"] + 1
      } else if (count >= 21 && count <= 30) {
        range_counts["21-30"] <- range_counts["21-30"] + 1
      } else if (count >= 31 && count <= 40) {
        range_counts["31-40"] <- range_counts["31-40"] + 1
      } else if (count >= 41 && count <= 50) {
        range_counts["41-50"] <- range_counts["41-50"] + 1
      } else if (count >= 51 && count <= 100) {
        range_counts["51-100"] <- range_counts["51-100"] + 1
      } else if (count >= 101 && count <= 150) {
        range_counts["101-150"] <- range_counts["101-150"] + 1
      } else if (count >= 151 && count <= 200) {
        range_counts["151-200"] <- range_counts["151-200"] + 1
      } else if (count >= 201 && count <= 1000) {
        range_counts["201-1000"] <- range_counts["201-1000"] + 1
      } else if (count > 1000) {
        range_counts[">1000"] <- range_counts[">1000"] + 1
      }
    }
    
    # Append the counts to the data frame
    for (range in names(range_counts)) {
      range_counts_df <<- rbind(range_counts_df, data.frame(Gene = gene_name, IdentifierType = type, Range = range, Count = range_counts[range]))
    }
  }
  
  # Count occurrences for canonical and non-canonical identifiers
  count_occurrences(canonical_identifiers, "canonical")
  count_occurrences(non_canonical_identifiers, "non-canonical")
}

# Write the output file
write.table(range_counts_df, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the range counts
print(range_counts_df)

# Print a message indicating completion
cat("Count ranges have been calculated and saved to", output_file, "\n")
