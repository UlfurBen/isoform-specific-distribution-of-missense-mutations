Step 1: Convert GTF to BED Format

# Load necessary libraries
library(dplyr)
library(stringr)

# Define output files
output_bed_file <- "Homo_sapiens.GRCh37.65.gtf"
output_bed_file_with_headers <- "Homo_sapiens.GRCh37.65_with_headers.bed"

# Read the GTF file
gtf_file <- "Homo_sapiens.GRCh37.65.bed"
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", quote="")

# Assign column names
colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Filter for "CDS" features
cds_data <- gtf_data %>% filter(feature == "CDS")

# Extract additional information from the attribute column
extract_attribute <- function(attribute, key) {
  str_match(attribute, paste0(key, " \"([^\"]+)\""))[, 2]
}

cds_data <- cds_data %>%
  mutate(gene_id = extract_attribute(attribute, "gene_id"),
         gene_name = extract_attribute(attribute, "gene_name"),
         name = gene_name)

# Select relevant columns for BED format
bed_df <- cds_data %>% select(chr, start, end, name, score, strand, frame, attribute)

# Write to a BED file without headers
write.table(bed_df, file=output_bed_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Add headers and write to a new file
colnames(bed_df) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "frame", "attribute")
write.table(bed_df, file=output_bed_file_with_headers, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("BED file with headers has been created.")


Step 2: Find Non-overlapping Regions

# Load necessary libraries
library(GenomicRanges)
library(dplyr)

# Define function to read BED file
read_bed_file <- function(file) {
  bed <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  bed$chromStart <- as.numeric(bed$chromStart)
  bed$chromEnd <- as.numeric(bed$chromEnd)
  bed <- bed %>% filter(!is.na(chromStart) & !is.na(chromEnd))
  return(bed)
}

# Define function to trim overlaps
trim_overlaps <- function(bed) {
  gr <- GRanges(seqnames = bed$chr,
                ranges = IRanges(start = bed$chromStart, end = bed$chromEnd),
                strand = bed$strand)
  
  overlaps <- findOverlaps(gr)
  for (i in seq_along(overlaps)) {
    query <- queryHits(overlaps[i])
    subject <- subjectHits(overlaps[i])
    
    if (query != subject) {
      q_start <- start(gr[query])
      q_end <- end(gr[query])
      s_start <- start(gr[subject])
      s_end <- end(gr[subject])
      
      if (q_start < s_start && q_end > s_start) end(gr[query]) <- s_start - 1
      if (s_start < q_start && s_end > q_start) end(gr[subject]) <- q_start - 1
    }
  }
  
  gr <- gr[start(gr) <= end(gr)]
  return(gr)
}

# Define function to write BED file
write_bed_file <- function(gr, output_file) {
  bed_out <- data.frame(seqnames = seqnames(gr),
                        chromStart = start(gr),
                        chromEnd = end(gr))
  write.table(bed_out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Process the BED file
input_file <- "Homo_sapiens.GRCh37.65_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.65_with_headers_isoform_specific_regions.bed"

bed <- read_bed_file(input_file)
gr_trimmed <- trim_overlaps(bed)
write_bed_file(gr_trimmed, output_file)

print("Trimmed BED file has been written.")


Step 3: Extract Specific Information and Filter

# Load necessary libraries
library(dplyr)
library(stringr)

# Define the file path
file_path <- "Homo_sapiens.GRCh37.65_with_headers_isoform_specific_regions.bed"

# Read the BED file into a data frame
bed_data <- read.table(file_path, header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t")

# Extract ENST and gene_name from the attribute column
extract_values <- function(attribute) {
  enst_value <- str_extract(attribute, "ENST\\S+")
  gene_name_value <- str_extract(attribute, "(?<=gene_name )\\S+")
  return(c(enst_value, gene_name_value))
}

# Apply the extraction function
values_extracted <- t(apply(bed_data, 1, function(row) extract_values(row["attribute"])))
bed_data <- cbind(bed_data, values_extracted)
colnames(bed_data)[(ncol(bed_data)-1):ncol(bed_data)] <- c("ENST", "gene_name")

# Filter to include only relevant columns and non-NA values
filtered_data <- bed_data %>%
  select(chr, chromStart, chromEnd, ENST, gene_name) %>%
  filter(!is.na(ENST) & !is.na(gene_name))

# Write to a new BED file
write.table(filtered_data, "filtered_Homo_sapiens.GRCh37.65.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

print("Filtering complete.")


Step 4: Sort and Filter BED File
                            
# Sort the BED file
sort -k1,1 -k2,2n filtered_Homo_sapiens.GRCh37.65.bed > sorted_Homo_sapiens.GRCh37.65_with_headers_isoform_specific_regions_bedtools_non_scientific.bed

Step 5: Remove Non-numeric and Sort Again

# Load necessary libraries
library(dplyr)

# Define function to check if a value is numeric or "MT", "X", "Y"
is_strictly_numeric_or_chr <- function(x) {
  grepl("^\\d+$", x) || x %in% c("MT", "X", "Y")
}

# Read and filter the BED file
input_file <- "sorted_Homo_sapiens.GRCh37.65_with_headers_isoform_specific_regions_bedtools_non_scientific.bed"
output_file <- "filtered_sorted_bed_file.bed"

bed_data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)
filtered_data <- bed_data %>%
  filter(sapply(V1, is_strictly_numeric_or_chr)) %>%
  mutate(V1 = factor(V1, levels = c(as.character(1:22), "X", "Y", "MT")),
         V2 = as.numeric(V2)) %>%
  arrange(V1, V2)

# Write to a new file
write.table(filtered_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("Filtering and sorting complete.")


Step 6: Remove Scientific Notation

# Load necessary libraries
library(dplyr)
library(readr)

# Define function to check for scientific notation
contains_scientific_notation <- function(row) {
  any(grepl("e\\+", row))
}

# Read and filter the BED file
file_path <- "filtered_sorted_bed_file.bed"
output_path <- "filtered_no_scientific_notation.bed"

bed_data <- read_delim(file_path, delim = "\t", col_names = FALSE)
filtered_data <- bed_data %>%
  filter(!apply(., 1, contains_scientific_notation))

# Write to a new file
write_delim(filtered_data, output_path, delim = "\t", col_names = FALSE)

print("Scientific notation filtering complete.")
