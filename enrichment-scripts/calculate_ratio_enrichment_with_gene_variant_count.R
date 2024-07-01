# Define file paths
input_file <- "merged_calculate_missense_variant_enrichment_with_gene_total_variant_count.txt"
output_file <- "merged_calculate_missense_variant_enrichment_ratio_with_gene_total_variant_count.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Convert necessary columns to numeric
data$Count <- as.numeric(data$Count)
data$Length <- as.numeric(data$Length)
data$Gene_Total_Count <- as.numeric(data$Gene_Total_Count)

# Calculate the ratio
data$Ratio <- (data$Count * data$Length) / data$Gene_Total_Count

# Write the modified data to the output file
write.table(data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the first few rows of the modified data to verify
print(head(data))

# Print a message indicating completion
cat("Data with calculated ratios has been saved to", output_file, "\n")
