# Load necessary libraries
library(dplyr)
library(data.table)

# Set the library path
# .libPaths("/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1")

# Define the input and output file paths
input_file_regions <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"
input_file_variants <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_with_variant_counts.bed"

# Read the input files
regions_data <- fread(input_file_regions, header = TRUE, sep = "\t")
variants_data <- fread(input_file_variants, header = TRUE, sep = "\t")

# Print column names for debugging
cat("Column names in regions_data:\n")
print(colnames(regions_data))
cat("\nColumn names in variants_data:\n")
print(colnames(variants_data))

# Convert 'chr' columns to character type to ensure consistency
regions_data[, chr := as.character(chr)]
variants_data[, chr := as.character(chr)]

# Create a key for fast joins
setkey(variants_data, chr, chromStart, chromEnd)

# Function to count variants in each region
count_variants <- function(chr, chromStart, chromEnd) {
  return(variants_data[chr == chr & chromStart <= chromEnd & chromEnd >= chromStart, .N])
}

# Apply the function to each row of the regions data
regions_data[, variant_count := mapply(count_variants, chr, chromStart, chromEnd)]

# Write the output to a new BED file with headers
fwrite(regions_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Variant counts have been added and data have been written to", output_file, "\n")
