# Load necessary libraries
library(dplyr)

# Define output files
output_bed_file <- "Homo_sapiens.GRCh37.87.bed"

# Read the GTF file
gtf_file <- "Homo_sapiens.GRCh37.87.gtf"

# Use system's grep command to filter lines with 'exon_id'
# This will return only lines containing 'exon_id'
exon_lines <- system(paste("grep 'exon_id' ", gtf_file), intern = TRUE)

# Write the grep output directly to the BED file
writeLines(exon_lines, con = output_bed_file)
