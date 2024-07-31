# Load necessary library
library(data.table)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Select only some columns
bed_data_subset <- bed_data[, .(V1, V2, V3, V5, V7, V9, V10 V16)]

# Write the modified data to a new file without headers
fwrite(bed_data_subset, output_file, sep = "\t", col.names = FALSE)

cat("File has been processed and saved as", output_file, "\n")





