# Read the data from the file
data <- read.table("homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract the second column
second_column <- data[, 2]

print(second_column)

# Calculate unique counts in the second column
unique_counts <- length(unique(second_column))

# Print the unique count
print(unique_counts)

