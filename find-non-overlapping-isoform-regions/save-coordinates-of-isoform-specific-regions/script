# Load necessary libraries
library(dplyr)
library(GenomicRanges)

# Define the input and output file paths
input_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff.bed"
output_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions.bed"

# Read the input file
data <- read.delim(input_file, header = TRUE, sep = "\t")

# Filter the data for entries with "CDS" in the score column
cds_data <- data %>% filter(score == "CDS")

# Replace non-standard strand values with '*'
cds_data$strand <- ifelse(cds_data$strand %in% c("+", "-", "*"), cds_data$strand, "*")

# Create a GRanges object
granges_data <- GRanges(seqnames = cds_data$chr,
                        ranges = IRanges(start = cds_data$chromStart,
                                         end = cds_data$chromEnd),
                        strand = cds_data$strand,
                        name = cds_data$name,
                        frame = cds_data$frame,
                        attribute = cds_data$attribute,
                        miscellaneous = cds_data$miscellaneous)

# Find non-overlapping regions
non_overlapping <- reduce(granges_data)

# Extract the data back into a data frame
non_overlapping_df <- data.frame(chr = seqnames(non_overlapping),
                                 chromStart = start(non_overlapping),
                                 chromEnd = end(non_overlapping),
                                 strand = strand(non_overlapping))

# Add back the original columns
non_overlapping_df <- non_overlapping_df %>%
  left_join(cds_data, by = c("chr", "chromStart", "chromEnd", "strand"))

# Select and order columns
final_output <- non_overlapping_df %>%
  select(chr, chromStart, chromEnd, name, score, strand, frame, attribute, miscellaneous)

# Write the output to a new BED file
write.table(final_output, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Non-overlapping regions have been written to", output_file, "\n")
