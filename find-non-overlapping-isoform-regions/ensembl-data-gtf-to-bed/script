# Load necessary library
library(dplyr)

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

# Swap columns: move 4th and 5th column to 2nd and 3rd positions and vice versa
bed_df <- cds_data %>% select(chr, start, end, source, feature, score, strand, frame, attribute)

# Write to a BED file
write.table(bed_df, file=output_bed_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Read the output BED file
bed_df <- read.table(output_bed_file, sep="\t", header=FALSE)

# Set column names
colnames(bed_df) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "frame", "attribute", "miscellaneous")

# Display the first 10 rows to confirm changes
print(head(bed_df, 10))

# Write the modified data frame back to a new file
write.table(bed_df, file=output_bed_file_with_headers, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
