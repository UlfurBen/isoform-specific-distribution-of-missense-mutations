# Set the library path
.libPaths("/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1")

# Load necessary libraries
library(dplyr)
library(data.table)

# Define the input and output file paths
input_file_regions <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions.bed"
input_file_variants <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions_with_variant_counts.bed"

# Read the input files
regions_data <- fread(input_file_regions, header = TRUE, sep = "\t")
variants_data <- fread(input_file_variants, header = TRUE, sep = "\t")

# Sort the variants data by chromosome and start position
variants_data <- variants_data[order(variants_data$chr, variants_data$chromStart)]

# Define a function to count the variants in each region
count_variants_in_region <- function(chr, start, end, variants_data) {
  subset_variants <- variants_data[chr == variants_data$chr & start <= variants_data$chromEnd & end >= variants_data$chromStart]
  return(nrow(subset_variants))
}

# Apply the function to each row of the regions data
regions_data <- regions_data %>%
  rowwise() %>%
  mutate(variant_count = count_variants_in_region(chr, chromStart, chromEnd, variants_data))

# Write the output to a new BED file
fwrite(regions_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Variant counts have been added and data have been written to", output_file, "\n")
