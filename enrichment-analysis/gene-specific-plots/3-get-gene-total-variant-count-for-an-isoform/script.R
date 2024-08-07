# Define file paths
input_file <- "enrichment_filtered_with_seq_len_more_than_before.txt"
output_file <- "enrichment_filtered_with_gene_total_variant_count.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Extract the gene identifier (part before the hyphen) from the Identifier column
data$Gene <- sub("-.*", "", data$Identifier)

# Calculate the sum of the Variation.count values for all isoforms of each gene
gene_counts <- aggregate(Variation.count ~ Gene, data = data, sum)

# Rename the Variation.count column to Gene_Total_Count
colnames(gene_counts)[2] <- "Gene_Total_Count"

# Merge the gene total counts back into the original data frame
merged_data <- merge(data, gene_counts, by = "Gene")

# Write the modified data to the output file
write.table(merged_data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the first few rows of the merged data to verify
print(head(merged_data))

# Print a message indicating completion
cat("Merged data with gene total counts has been saved to", output_file, "\n")
