# Define file paths
input_file <- "calculate_ratio_with_isoform_gene_count.txt"
output_file <- "count_ranges.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Convert Count column to numeric, handling potential non-numeric entries
data$Count <- as.numeric(data$Count)

# Define the count ranges and initialize counts
ranges <- c("0-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-150", "151-200", "201-1000", ">1000")
range_counts <- setNames(rep(0, length(ranges)), ranges)

# Count the number of lines that fall into each range
range_counts["0-5"] <- sum(data$Count >= 0 & data$Count <= 5, na.rm = TRUE)
range_counts["6-10"] <- sum(data$Count >= 6 & data$Count <= 10, na.rm = TRUE)
range_counts["11-20"] <- sum(data$Count >= 11 & data$Count <= 20, na.rm = TRUE)
range_counts["21-30"] <- sum(data$Count >= 21 & data$Count <= 30, na.rm = TRUE)
range_counts["31-40"] <- sum(data$Count >= 31 & data$Count <= 40, na.rm = TRUE)
range_counts["41-50"] <- sum(data$Count >= 41 & data$Count <= 50, na.rm = TRUE)
range_counts["51-100"] <- sum(data$Count >= 51 & data$Count <= 100, na.rm = TRUE)
range_counts["101-150"] <- sum(data$Count >= 101 & data$Count <= 150, na.rm = TRUE)
range_counts["151-200"] <- sum(data$Count >= 151 & data$Count <= 200, na.rm = TRUE)
range_counts["201-1000"] <- sum(data$Count >= 201 & data$Count <= 1000, na.rm = TRUE)
range_counts[">1000"] <- sum(data$Count > 1000, na.rm = TRUE)

# Convert the range_counts to a data frame
range_counts_df <- data.frame(Range = names(range_counts), Count = unlist(range_counts))

# Write the output file
write.table(range_counts_df, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the range counts
print(range_counts_df)

# Print a message indicating completion
cat("Count ranges have been calculated and saved to", output_file, "\n")

