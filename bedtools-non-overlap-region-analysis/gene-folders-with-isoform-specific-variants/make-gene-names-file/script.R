# Extract the Gene_Name column
gene_names <- data$Gene_Name

# Write the gene names to a CSV file without quotation marks
write.csv(gene_names, file = "gene_names.csv", row.names = FALSE, quote = FALSE)
