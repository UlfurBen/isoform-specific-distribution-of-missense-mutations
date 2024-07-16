# Define file paths
input_file <- "enrichment_filtered_with_gene_total_variant_count.txt"
output_file <- "count_ranges_with_zero.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Check column names and correct if necessary
colnames(data) <- trimws(colnames(data))  # Remove any leading/trailing spaces from column names
expected_colname <- "Variation.count"
if (!expected_colname %in% colnames(data)) {
  stop(paste("Expected column", expected_colname, "not found in input file."))
}

# Convert Variation.count column to numeric, handling potential non-numeric entries
data$Variation.count <- as.numeric(data$Variation.count)

# Define the count ranges and initialize counts, including a count for zero variants
ranges <- c("0", "0-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-150", "151-200", "201-1000", ">1000")
range_counts <- setNames(rep(0, length(ranges)), ranges)

# Count the number of lines that fall into each range
range_counts["0"] <- sum(data$Variation.count == 0, na.rm = TRUE)
range_counts["0-5"] <- sum(data$Variation.count >= 1 & data$Variation.count <= 5, na.rm = TRUE)
range_counts["6-10"] <- sum(data$Variation.count >= 6 & data$Variation.count <= 10, na.rm = TRUE)
range_counts["11-20"] <- sum(data$Variation.count >= 11 & data$Variation.count <= 20, na.rm = TRUE)
range_counts["21-30"] <- sum(data$Variation.count >= 21 & data$Variation.count <= 30, na.rm = TRUE)
range_counts["31-40"] <- sum(data$Variation.count >= 31 & data$Variation.count <= 40, na.rm = TRUE)
range_counts["41-50"] <- sum(data$Variation.count >= 41 & data$Variation.count <= 50, na.rm = TRUE)
range_counts["51-100"] <- sum(data$Variation.count >= 51 & data$Variation.count <= 100, na.rm = TRUE)
range_counts["101-150"] <- sum(data$Variation.count >= 101 & data$Variation.count <= 150, na.rm = TRUE)
range_counts["151-200"] <- sum(data$Variation.count >= 151 & data$Variation.count <= 200, na.rm = TRUE)
range_counts["201-1000"] <- sum(data$Variation.count >= 201 & data$Variation.count <= 1000, na.rm = TRUE)
range_counts[">1000"] <- sum(data$Variation.count > 1000, na.rm = TRUE)

# Convert the range_counts to a data frame
range_counts_df <- data.frame(Range = names(range_counts), Count = unlist(range_counts))

# Write the output file
write.table(range_counts_df, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the range counts
print(range_counts_df)

# Print a message indicating completion
cat("Count ranges have been calculated and saved to", output_file, "\n")
