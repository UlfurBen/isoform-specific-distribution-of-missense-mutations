# Load necessary library
# Install the library if not already installed
if (!require(readr)) install.packages("readr", dependencies=TRUE)

# Read the BED file (assuming tab-delimited format)
bed_file <- readr::read_tsv("variants_benign_pathogenic_non_vus_non_conflicting_including_non_annotation.bed", col_names = FALSE)

# Extract unique values from the 5th column
unique_values <- unique(bed_file[[5]])

# Write the unique values to a CSV file
write.csv(unique_values, "whole_genome.csv", row.names = FALSE, quote = FALSE)

print("Unique values saved to whole_genome.csv")
