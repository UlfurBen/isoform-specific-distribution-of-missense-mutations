# Load required libraries
library(utils)

# Define file paths
# Homo_sapiens_variation.txt.gz is on elja (too big for github)
input_file <- "homo_sapiens_variation.txt.gz"
output_file <- "homo_sapiens_variation_missense_ClinVar.txt"

# Read the input file
data <- read.table(gzfile(input_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter the data to only include lines that contain both 'missense' and 'ClinVar'
filtered_data <- data[grepl("missense", data$V1) & grepl("ClinVar", data$V1), ]

# Write the filtered data to the output file
write.table(filtered_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print a message indicating completion
cat("Filtered data has been saved to", output_file, "\n")
