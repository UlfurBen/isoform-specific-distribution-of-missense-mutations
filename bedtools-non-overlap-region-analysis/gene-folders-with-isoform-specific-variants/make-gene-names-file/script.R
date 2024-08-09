# Load the data from The-Epigenetic-Machinery.csv file with headers
data <- read.csv("The-Epigenetic-Machinery.csv", stringsAsFactors = FALSE)

# Extract the Gene_Name column
gene_names <- data$Gene_Name

# Write the gene names to a CSV file without quotation marks and with headers
write.csv(gene_names, file = "gene_names.csv", row.names = FALSE, quote = FALSE)

# Manually remove first column from gene_names.csv
