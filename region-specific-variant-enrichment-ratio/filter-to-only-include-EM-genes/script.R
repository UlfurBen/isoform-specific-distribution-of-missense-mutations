data.table_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/data.table"
dplyr_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/dplyr"

# Add these paths to the library paths
.libPaths(c(dirname(data.table_path), .libPaths()))
.libPaths(c(dirname(dplyr_path), .libPaths()))

# Load necessary libraries
library(data.table)
library(dplyr)

# Define file paths for input and output files
epigenetic_file <- "The-Epigenetic-Machinery.csv"
bed_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_with_variant_counts_enrichment_ratio_in_descending_order.bed"
output_file <- "Matched_Genes.bed"

# Read the CSV file into a data table
epigenetic_data <- fread(epigenetic_file)

# Read the BED file into a data table
bed_data <- fread(bed_file, header = TRUE)

# Extract the Gene_Name column from the epigenetic data
gene_names <- epigenetic_data$Gene_Name

# Initialize an empty data table to store matched rows
matched_genes <- data.table()

# Iterate over each gene name and search the entire row in the BED file
for (gene in gene_names) {
  # Filter rows in bed_data that contain the gene name in any part of the row
  matched_rows <- bed_data[apply(bed_data, 1, function(row) any(grepl(gene, row)))]
  
  # Append matched rows to matched_genes
  matched_genes <- rbind(matched_genes, matched_rows)
}

# Remove duplicate rows (if any)
matched_genes <- unique(matched_genes)

# Save the matched rows to a new file
fwrite(matched_genes, file = output_file, sep = "\t", col.names = TRUE)

# Display the column names of the matched genes data
print(colnames(matched_genes))
