# Load necessary library
library(dplyr)
library(stringr)

# Define output files
output_bed_file <- "Homo_sapiens.GRCh37.87.bed"
output_bed_file_with_headers <- "Homo_sapiens.GRCh37.87_with_headers.bed"

# Read the GTF file
gtf_file <- "Homo_sapiens.GRCh37.87.gtf"
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", quote="")

# Assign column names, initializing the first column as "chr"
colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Filter the data to include only lines containing "CDS"
cds_data <- gtf_data %>% filter(feature == "CDS")

# Extract additional information from the attribute column
extract_attribute <- function(attribute, key) {
  str_match(attribute, paste0(key, " \"([^\"]+)\""))[, 2]
}

cds_data <- cds_data %>%
  mutate(gene_id = extract_attribute(attribute, "gene_id"),
         gene_version = extract_attribute(attribute, "gene_version"),
         gene_name = extract_attribute(attribute, "gene_name"),
         gene_source = extract_attribute(attribute, "gene_source"),
         gene_biotype = extract_attribute(attribute, "gene_biotype"),
         name = gene_name,  # For simplicity, using gene_name as 'name'
         miscellaneous = NA)  # Placeholder for miscellaneous

# Rearrange columns for the BED file format
bed_df <- cds_data %>% select(chr, start, end, name, score, strand, frame, attribute, miscellaneous)

# Write to a BED file
write.table(bed_df, file=output_bed_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Read the output BED file
bed_df <- read.table(output_bed_file, sep="\t", header=FALSE)

# Set column names correctly
colnames(bed_df) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "frame", "attribute", "miscellaneous")

# Display the first 10 rows to confirm changes
print(head(bed_df, 10))

# Write the modified data frame back to a new file with headers
write.table(bed_df, file=output_bed_file_with_headers, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
