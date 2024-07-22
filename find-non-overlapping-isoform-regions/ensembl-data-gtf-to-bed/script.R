# Load necessary library
library(dplyr)
library(tidyr)

# Define output files
output_bed_file <- "Homo_sapiens.GRCh37.87.bed"
output_bed_file_with_headers <- "Homo_sapiens.GRCh37.87_with_headers.bed"

# Read the GTF file, skipping lines starting with '#!'
gtf_file <- "Homo_sapiens.GRCh37.87.gtf"
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", quote="")

# Assign column names, initializing the first column as "chr"
colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Filter the data to include only lines containing "CDS"
cds_data <- gtf_data %>% filter(feature == "CDS")

# Add additional columns with NA values
cds_data <- cds_data %>%
  mutate(
    gene_id = NA, gene_version = NA, transcript_id = NA, transcript_version = NA,
    exon_number = NA, gene_name = NA, gene_source = NA, gene_biotype = NA,
    transcript_name = NA, transcript_source = NA, transcript_biotype = NA,
    protein_id = NA, protein_version = NA
  )

# Select relevant columns and rearrange
bed_df <- cds_data %>%
  select(chr, start, end, source, feature, score, strand, frame, attribute, 
         gene_id, gene_version, transcript_id, transcript_version, exon_number, 
         gene_name, gene_source, gene_biotype, transcript_name, transcript_source, 
         transcript_biotype, protein_id, protein_version)

# Write to a BED file
write.table(bed_df, file=output_bed_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Read the output BED file
bed_df <- read.table(output_bed_file, sep="\t", header=FALSE)

# Set column names
colnames(bed_df) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "frame", "attribute", 
                      "gene_id", "gene_version", "transcript_id", "transcript_version", "exon_number", 
                      "gene_name", "gene_source", "gene_biotype", "transcript_name", "transcript_source", 
                      "transcript_biotype", "protein_id", "protein_version")

# Display the first 10 rows to confirm changes
print(head(bed_df, 10))

# Write the modified data frame back to a new file with headers
write.table(bed_df, file=output_bed_file_with_headers, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
