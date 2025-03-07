# Specify the data.table path on elja
data.table_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/data.table"

# Add data.table_path to the library path
.libPaths(c(dirname(data.table_path), .libPaths()))

# Load necessary libraries
library(data.table)

# Define the input and output file paths
input_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_with_variant_counts_enrichment_ratio.bed"
output_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_with_variant_counts_enrichment_ratio_in_descending_order.bed"

# Read the input file
data <- fread(input_file, header = FALSE, sep = "\t")

# Print the column names for debugging
cat("Column names in the input file:\n")
print(colnames(data))

# Exclude rows where the second and third column values are the same
data <- data[V2 != V3]

# Order the rows in descending order based on the last column values
data <- data[order(-V11)]

# Write the output to a new file
fwrite(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Data has been ordered and written to", output_file, "\n")
