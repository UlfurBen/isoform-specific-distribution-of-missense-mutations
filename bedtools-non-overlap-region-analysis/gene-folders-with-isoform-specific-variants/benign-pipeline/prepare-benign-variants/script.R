# This script prepares the variant file for the benign isoform specific variant count analysis

#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(dplyr)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file_1 <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_benign.bed"
output_file_2 <- "variants-sorted-benign.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Define the exclusion terms
exclusion_terms <- c("pathogenic", "variant of uncertain significance", "conflicting")

# Filter the data to include only rows where 'Benign' appears in column 10 (V10), case insensitive,
# and exclude rows with conflicting terms in the same column
filtered_data <- bed_data %>%
  filter(grepl("benign", V10, ignore.case = TRUE) & !grepl(paste(exclusion_terms, collapse = "|"), V10, ignore.case = TRUE))

# Write the filtered data to a new file
fwrite(filtered_data, output_file_1, sep = "\t", col.names = FALSE)




# This next part of the script processes the benign isoform specific variant data and creates a sorted output file

# Read the benign BED file
bed_data_2 <- fread(output_file_1, header = FALSE)

# Extract the string that starts with 'ENST' from each row and add it as a new column (V6)
bed_data_2 <- bed_data_2 %>%
  mutate(V6 = sapply(V16, function(x) {
    matches <- regmatches(x, regexpr("ENST\\w+", x))
    ifelse(length(matches) > 0, matches, NA)
  }))

# Select the first 5 columns and the new 6th column
output_data <- bed_data_2 %>%
  select(V1:V5, V6)

# Write the processed data to a new file
fwrite(output_data, output_file_2, sep = "\t", col.names = FALSE)
cat("Sorted file with ENST strings has been processed and saved as", output_file_2, "\n")
