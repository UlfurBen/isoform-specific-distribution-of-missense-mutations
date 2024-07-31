# Load necessary library
library(data.table)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Select only some columns
bed_data_subset <- bed_data[, .(V1, V2, V3, V5, V7, V9, V10, V16)]

# Write the modified data to a new file without headers
fwrite(bed_data_subset, output_file, sep = "\t", col.names = FALSE)

cat("File has been processed and saved as", output_file, "\n")

# Load the necessary library
library(dplyr)

# Read the BED file into a data frame
# Replace 'your_file_path' with the actual path to your BED file
bed_data <- read.delim('homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed', header = FALSE)

# Define the exclusion terms
exclusion_terms <- c("variant of uncertain significance", "benign")

# Filter the data to include only 'Pathogenic' in column 7 and exclude the specified terms
filtered_data <- bed_data %>%
  filter(tolower(V7) == "pathogenic" & !tolower(V7) %in% exclusion_terms)

# Write the filtered data to a new file
write.table(filtered_data, 'homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_pathogenic.bed', sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)






