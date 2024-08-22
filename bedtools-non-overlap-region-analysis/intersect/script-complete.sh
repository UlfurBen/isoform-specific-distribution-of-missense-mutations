# Filter input bed file to include first 3 columns and two values in column 9

# Load necessary library
library(dplyr)
library(stringr)

# Define the file path (adjust as necessary)
file_path <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"

# Read the BED file into a data frame
# Specify only the first 9 columns, assuming there are no headers
bed_data <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", colClasses = "character")

# Function to extract ENST and gene_name values from column 9
extract_values <- function(col9) {
  enst_value <- str_extract(col9, "ENST\\S+")
  gene_name_value <- str_extract(col9, "(?<=gene_name )\\S+")
  # Remove trailing semicolon from both values if they exist
  enst_value <- str_replace(enst_value, ";$", "")
  gene_name_value <- str_replace(gene_name_value, ";$", "")
  return(c(enst_value, gene_name_value))
}

# Apply the function to column 9 and create new columns
values_extracted <- t(apply(bed_data, 1, function(row) extract_values(row[9])))
bed_data <- cbind(bed_data, values_extracted)
colnames(bed_data)[(ncol(bed_data)-1):ncol(bed_data)] <- c("ENST", "gene_name")

# Filter the data to include only the first 3 columns and the extracted ENST and gene_name columns
filtered_data <- bed_data %>%
  select(V1, V2, V3, ENST, gene_name) %>%
  filter(!is.na(ENST) & !is.na(gene_name))

# Write the filtered data to a new BED file
write.table(filtered_data, "filtered_Homo_sapiens.GRCh37.87.bed", 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Print a message indicating the script has finished
cat("Filtering complete. The filtered file has been saved as 'filtered_Homo_sapiens.GRCh37.87.bed'.\n")







# Correct input files

# filtered_Homo_sapiens.GRCh37.87.bed


# Sort the variant bed file
sort -k1,1 -k2,2n homo_sapiens_variation_missense_ClinVar_filtered_relevancy_pathogenic.bed > sorted_hosmo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_pathogenic.bed

# Sort the region bed file
sort -k1,1 -k2,2n filtered_Homo_sapiens.GRCh37.87.bed > sorted_Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_bedtools_non_scientific.bed






# Remove all rows with first column value that isn't numeric and sort the file

# Load necessary library
library(dplyr)

# Function to check if a value is numeric or "MT", "X", "Y" without coercion
is_strictly_numeric_or_chr <- function(x) {
  grepl("^\\d+$", x) || x %in% c("MT", "X", "Y")
}

# Read the BED file
bed_file <- "sorted_Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_bedtools_non_scientific.bed"
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)

# Filter rows where the first column is strictly numeric or "MT", "X", "Y"
filtered_data <- bed_data %>%
  filter(sapply(V1, is_strictly_numeric_or_chr))

# Convert first two columns to appropriate types for sorting
filtered_data <- filtered_data %>%
  mutate(V1 = factor(V1, levels = c(as.character(1:22), "X", "Y", "MT")), V2 = as.numeric(V2))

# Sort the data by the first two columns
sorted_data <- filtered_data %>%
  arrange(V1, V2)

# Write the sorted data back to a new BED file
write.table(sorted_data, "filtered_sorted_bed_file.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print("Filtering and sorting complete. The new file is saved as 'filtered_sorted_bed_file.bed'.")







# Remove scientific notation

# Load necessary library
library(dplyr)
library(readr)

# Define the file path (adjust as necessary)
file_path <- "filtered_sorted_bed_file.bed"
output_path <- "filtered_no_scientific_notation.bed"

# Read the BED file into a data frame
bed_data <- read_delim(file_path, delim = "\t", col_names = FALSE)

# Function to check for scientific notation in any column
contains_scientific_notation <- function(row) {
  any(grepl("e\\+", row))
}

# Filter out rows that contain scientific notation
filtered_data <- bed_data %>% 
  filter(!apply(., 1, contains_scientific_notation))

# Write the filtered data to a new BED file
write_delim(filtered_data, output_path, delim = "\t", col_names = FALSE)

# Print a message indicating the script has finished
cat("Filtering complete. The filtered file has been saved as 'filtered_no_scientific_notation.bed'.\n")




# Find non-overlap regions

bedtools intersect -a filtered_no_scientific_notation.bed -b filtered_no_scientific_notation.bed -wo | awk '$2 != $7 || $3 != $8 || $4 != $9' > partial_overlapping.bed

bedtools subtract -a filtered_no_scientific_notation.bed -b partial_overlapping.bed > shortened_regions.bed


# Find variants count in each region

# shortened_regions.bed
# sorted_hosmo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_pathogenic.bed

bedtools intersect -a shortened_regions.bed -b sorted_hosmo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_pathogenic.bed -c > intersected_variants.bed






# Sort rows by variant count value from highest to lowest

sort -k6,6nr intersected_variants.bed > intersected_variants_sorted.bed

# Format

column -t intersected_variants_sorted.bed > intersected_variants_sorted_formatted_lp_and_p.bed



# Calculate enrichment ratio

awk '{diff = $3 - $2; enrichment_ratio = (diff == 0) ? 0 : $6 / diff; print $0, enrichment_ratio}' intersected_variants_sorted.bed > intersected_variants_with_enrichment.bed


sort -k7,7nr intersected_variants_with_enrichment.bed > intersected_variants_with_enrichment_sorted_lp_p.bed



# Filter to only include enrichment ratio >0

# Likely pathogenic and pathogenic variants
awk '$7 != 0' intersected_variants_with_enrichment_sorted_lp_p.bed | sort -k7,7nr > intersected_variants_filtered_sorted_with_X_Y_MT_pathogenic_and_likely_pathogenic_variant_count.bed
column -t intersected_variants_filtered_sorted_with_X_Y_MT_pathogenic_and_likely_pathogenic_variant_count.bed > formatted_variants_likely_p_p.bed






----------------------------------------------------------------------------------------------------------------------





# Only pathogenic variants
awk '$7 != 0' intersected_variants_with_enrichment_sorted.bed | sort -k7,7nr > intersected_variants_filtered_sorted_with_X_Y_MT_pathogenic_variant_count.bed
column -t intersected_variants_filtered_sorted_with_X_Y_MT_pathogenic_variant_count.bed > formatted_variants.bed

# New command
awk '$7 != 0' intersected_variants_with_enrichment_sorted.bed | sort -k7,7nr > intersected_variants_filtered_sorted_with_X_Y_MT.bed
# Old command
awk '$7 != 0' intersected_variants_with_enrichment_sorted.bed | sort -k7,7nr > intersected_variants_filtered_sorted.bed
