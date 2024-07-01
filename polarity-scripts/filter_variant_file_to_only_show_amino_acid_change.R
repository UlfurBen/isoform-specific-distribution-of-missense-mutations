# Define file paths
input_file <- "filter_database_by_EM_genes_missense_variants_ClinVar_with_RCV.txt"
output_file <- "filter_input_to_only_contain_categorical_change_for_each_variant.txt"

# Read the input file
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# Function to remove numeric characters between amino acid names
filter_amino_acids <- function(variant) {
  gsub("p\\.([A-Za-z]+)\\d+([A-Za-z]+)", "\\1\\2", variant)
}

# Apply the function to the relevant column (assuming the third column contains the amino acid changes)
data[, 3] <- sapply(data[, 3], filter_amino_acids)

# Write the modified data to the output file
write.table(data, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Print a message indicating completion
cat("Filtered data has been saved to", output_file, "\n")
