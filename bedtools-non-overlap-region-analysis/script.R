data.table_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/data.table"
stringr_path <- "/hpcapps/lib-mimir/software/R/4.1.2-foss-2021b/lib64/R/library/stringr"

.libPaths(c(dirname(data.table_path), .libPaths()))
.libPaths(c(dirname(stringr_path), .libPaths()))

# Load necessary libraries
library(data.table) # For data manipulation
library(stringr)    # For string operations

# Define file paths
region_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
variant_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"

# Read the region and variant files
regions <- fread(region_file, sep = "\t", header = TRUE)
variants <- fread(variant_file, sep = "\t", header = TRUE)

# Save regions and variants to temporary files for bedtools usage
write.table(regions, "regions.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(variants, "variants.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Use bedtools to merge overlapping regions and split into non-overlapping regions
system("bedtools merge -i regions.bed > merged_regions.bed")
system("bedtools complement -i merged_regions.bed -g <(cut -f1,2 regions.bed | sort | uniq) > non_overlapping_regions.bed")

# Read the non-overlapping regions back into R
non_overlapping_regions <- fread("non_overlapping_regions.bed", sep = "\t", header = FALSE)

# Rename columns for clarity
setnames(non_overlapping_regions, c("chr", "chromStart", "chromEnd"))

# Use bedtools to count variants in non-overlapping regions
system("bedtools intersect -a non_overlapping_regions.bed -b variants.bed -c > variant_counts.bed")

# Read the variant counts back into R
variant_counts <- fread("variant_counts.bed", sep = "\t", header = FALSE)

# Rename columns for clarity
setnames(variant_counts, c("chr", "chromStart", "chromEnd", "variant_count"))

# Display the result
print(variant_counts)

# Clean up temporary files
file.remove("regions.bed", "variants.bed", "merged_regions.bed", "non_overlapping_regions.bed", "variant_counts.bed")
