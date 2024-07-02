# Define file paths
input_file <- "merged_calculate_missense_variant_enrichment_within_isoforms_with_lengths_all_300.txt"
output_file <- "merged_calculate_missense_variant_enrichment_with_gene_total_variant_count.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Extract the gene identifier (part before the hyphen) from the Isoform column
data$Gene <- sub("-.*", "", data$Isoform)

# Calculate the sum of the Count values for all isoforms of each gene
gene_counts <- aggregate(Count ~ Gene, data = data, sum)

# Rename the Count column to Gene_Total_Count
colnames(gene_counts)[2] <- "Gene_Total_Count"

# Merge the gene total counts back into the original data frame
merged_data <- merge(data, gene_counts, by = "Gene")

# Write the modified data to the output file
write.table(merged_data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the first few rows of the merged data to verify
print(head(merged_data))

# Print a message indicating completion
cat("Merged data with gene total counts has been saved to", output_file, "\n")
