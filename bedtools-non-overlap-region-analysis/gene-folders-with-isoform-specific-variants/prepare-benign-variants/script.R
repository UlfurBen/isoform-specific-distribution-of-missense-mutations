# This script prepares the variant file for the benign isoform specific variant count analysis

#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(dplyr)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_benign.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Define the exclusion terms
exclusion_terms <- c("pathogenic", "variant of uncertain significance", "conflicting")

# Filter the data to include only rows where 'Benign' appears in column 10 (V10), case insensitive,
# and exclude rows with conflicting terms in the same column
filtered_data <- bed_data %>%
  filter(grepl("benign", V10, ignore.case = TRUE) & !grepl(paste(exclusion_terms, collapse = "|"), V10, ignore.case = TRUE))

# Write the filtered data to a new file
fwrite(filtered_data, output_file, sep = "\t", col.names = FALSE)
cat("Filtered file has been processed and saved as", output_file, "\n")
