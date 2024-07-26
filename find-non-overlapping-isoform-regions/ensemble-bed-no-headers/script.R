# Load necessary library
library(data.table)

# Define file paths
input_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.87_without_headers.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Select only the first three columns (chr, chromStart, chromEnd)
bed_data_subset <- bed_data[, .(V1, V2, V3)]

# Filter out contigs and unplaced sequences
standard_chromosomes <- c(as.character(1:22), "X", "Y", "MT")
bed_data_filtered <- bed_data_subset[V1 %in% standard_chromosomes]

# Write the modified data to a new file without headers
fwrite(bed_data_filtered, output_file, sep = "\t", col.names = FALSE)

cat("File has been processed and saved as", output_file, "\n")
