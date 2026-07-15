# ==============================================================================
# Purpose: Identify genomic regions (exons) with high pathogenic variant density.
# Method: Merges genomic regions, computes non-overlapping gaps, and counts variants.
# ==============================================================================

# 1. SETUP & DEPENDENCIES -----------------------------------------------------
library(data.table)

# Prevent scientific notation in genomic coordinates
options(scipen = 999)

# Define descriptive file paths
FILE_EXON_REGIONS  <- "Homo_sapiens.GRCh37.87_with_headers.bed"
FILE_CLINVAR_VARS  <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
FILE_FINAL_OUTPUT  <- "exon_pathogenic_variant_density.tsv"

# 2. LOAD & PREPARE DATA -------------------------------------------------------
# Read datasets using data.table's fast fread
regions  <- fread(FILE_EXON_REGIONS, sep = "\t", header = TRUE)
variants <- fread(FILE_CLINVAR_VARS, sep = "\t", header = TRUE)

# Ensure data is sorted by chromosome and starting position
setorder(regions, chr, chromStart)
setorder(variants, chr, chromStart)

# 3. MERGE OVERLAPPING REGIONS (Replaces 'bedtools merge') ---------------------
# Reduce overlapping exons into unified, contiguous genomic blocks
merged_regions <- regions |>
  _[, .(chromStart = min(chromStart), chromEnd = max(chromEnd)), 
    by = .(chr, cumsum(chromStart > shift(chromEnd, fill = FIRST(chromStart))))]

# 4. COMPUTE NON-OVERLAPPING REGIONS (Replaces 'bedtools complement') ---------
# Identify the gaps between the merged exon regions
non_overlapping <- merged_regions |>
  _[, .(
    chromStart = chromEnd[-.N], 
    chromEnd   = chromStart[-1]
  ), by = chr] |>
  _[chromStart < chromEnd] # Filter out invalid or zero-length gaps

# 5. COUNT VARIANTS IN GAPS (Replaces 'bedtools intersect -c') -----------------
# Set keys to perform a highly optimized data.table overlap join
setkey(variants, chr, chromStart, chromEnd)

# Count how many variants fall strictly within each non-overlapping gap
variant_counts <- foverlaps(non_overlapping, variants, type = "within", nomatch = NULL) |>
  _[, .(variant_count = .N), by = .(chr, chromStart = i.chromStart, chromEnd = i.chromEnd)]

# Bring back any non-overlapping regions that had zero variants
variant_counts <- merge(non_overlapping, variant_counts, by = c("chr", "chromStart", "chromEnd"), all.x = TRUE)
variant_counts[is.na(variant_count), variant_count := 0]

# 6. CALCULATE DENSITY & SORT --------------------------------------------------
# Compute variants per base pair, safely handling potential division by zero
variant_counts[, region_length := chromEnd - chromStart]
variant_counts[, variant_density := variant_count / fcase(region_length == 0, 1, default = region_length)]

# Rank regions by highest mutation density first
setorder(variant_counts, -variant_density)

# 7. SAVE RESULTS --------------------------------------------------------------
# Export clean, tab-separated results
fwrite(variant_counts[, .(chr, chromStart, chromEnd, variant_count, variant_density)], 
       file = FILE_FINAL_OUTPUT, sep = "\t")

cat("Process complete. Output saved to:", FILE_FINAL_OUTPUT, "\n")
