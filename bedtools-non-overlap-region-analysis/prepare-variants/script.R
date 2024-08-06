#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(dplyr)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_pathogenic.bed"
final_output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_pathogenic_fresh.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Define the exclusion terms
exclusion_terms <- c("benign", "variant of uncertain significance")

# Filter the data to include only rows where 'Pathogenic' appears in column 10 (V10), case insensitive,
# and exclude rows with conflicting terms in the same column
filtered_data <- bed_data %>%
  filter(grepl("pathogenic", V10, ignore.case = TRUE) & !grepl(paste(exclusion_terms, collapse = "|"), V10, ignore.case = TRUE))

# Write the filtered data to a new file
fwrite(filtered_data, output_file, sep = "\t", col.names = FALSE)
cat("Filtered file has been processed and saved as", output_file, "\n")

# Extract specific columns: 1st, 2nd, 3rd, 4th, 5th, 10th, and 16th
selected_columns <- filtered_data %>% select(V1, V2, V3, V4, V5, V10, V16)

# Write the selected columns to the final output file
write.table(selected_columns, final_output_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
cat("Selected columns have been saved as", final_output_file, "\n")
