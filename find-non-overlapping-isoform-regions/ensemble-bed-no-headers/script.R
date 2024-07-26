# Load necessary library
library(data.table)

# Read the BED file
input_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.87_without_headers.bed"

# Read the file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Select only the first three columns (chr, chromStart, chromEnd)
bed_data_subset <- bed_data[, .(V1, V2, V3)]

# Write the modified data to a new file without headers
fwrite(bed_data_subset, output_file, sep = "\t", col.names = FALSE)

cat("File has been processed and saved as", output_file, "\n")
