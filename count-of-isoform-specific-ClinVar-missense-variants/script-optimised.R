# Load necessary libraries
library(dplyr)
library(data.table)

# Set the library path
.libPaths("/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1")

# Define the input and output file paths
input_file_regions <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions.bed"
input_file_variants <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions_with_variant_counts.bed"

# Read the input files
regions_data <- fread(input_file_regions, header = TRUE, sep = "\t")
variants_data <- fread(input_file_variants, header = TRUE, sep = "\t")

# Rename columns for merging
setnames(variants_data, old = c("chr", "chromStart", "chromEnd"), new = c("chr", "start", "end"))

# Create a key for fast joins
setkey(variants_data, chr, start, end)

# Function to count variants in each region
count_variants <- function(chr, start, end) {
  variants_in_region <- variants_data[chr == chr & start <= end & end >= start, .N]
  return(variants_in_region)
}

# Apply the function to each row of the regions data
regions_data[, variant_count := count_variants(chr, chromStart, chromEnd), by = .(chr, chromStart, chromEnd)]

# Write the output to a new BED file
fwrite(regions_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Variant counts have been added and data have been written to", output_file, "\n")
