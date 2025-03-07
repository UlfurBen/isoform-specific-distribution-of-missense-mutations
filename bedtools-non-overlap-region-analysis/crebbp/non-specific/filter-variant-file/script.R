# Filter input bed file to include first 3 columns and two values in column 9

# Load necessary library
library(dplyr)
library(stringr)

# Define the file path (adjust as necessary)
file_path <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy_pathogenic.bed"

# Read the BED file into a data frame
# Specify only the first 9 columns, assuming there are no headers
bed_data <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", colClasses = "character")

# Filter the data to include only the first 3 columns and the extracted ENST and gene_name columns
filtered_data <- bed_data %>%
  select(V1, V2, V3, V4, V5, V6, V16)

# Write the filtered data to a new BED file
write.table(filtered_data, "variants_filtered.bed", 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Print a message indicating the script has finished
cat("Filtering complete. The filtered file has been saved as 'variants_filtered.bed'.\n")

# Sort the variant bed file on command line
# sort -k1,1 -k2,2n variants_filtered.bed > variants_filtered_sorted.bed
