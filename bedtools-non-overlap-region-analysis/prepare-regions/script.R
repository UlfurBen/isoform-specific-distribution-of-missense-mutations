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





input_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"

# Set the library paths
# library_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1"
# .libPaths(library_path)

GenomicRanges_path <- "/hpcapps/lib-mimir/software/R/4.1.2-foss-2021b/lib64/R/library/GenomicRanges"
dplyr_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/dplyr"

# Print the paths to verify
print(GenomicRanges_path)
print(dplyr_path)

# Add these paths to the library paths
.libPaths(c(dirname(GenomicRanges_path), .libPaths()))
.libPaths(c(dirname(dplyr_path), .libPaths()))

# Load the packages
library(dplyr)
library(GenomicRanges)

# Define functions
read_bed_file <- function(file) {
  bed <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end", "source", "feature", "score", "strand", "frame", "attribute")
  
  bed$start <- suppressWarnings(as.numeric(bed$start))
  bed$end <- suppressWarnings(as.numeric(bed$end))
  bed <- bed %>% filter(!is.na(start) & !is.na(end))
  
  return(list(data = bed, headers = colnames(bed)))
}

trim_overlaps <- function(bed) {
  gr <- GRanges(seqnames = bed$chr,
                ranges = IRanges(start = bed$start, end = bed$end),
                strand = bed$strand,
                source = bed$source,
                feature = bed$feature,
                score = bed$score,
                frame = bed$frame,
                attribute = bed$attribute)
  
  overlaps <- findOverlaps(gr)
  
  for (i in seq_along(overlaps)) {
    query <- queryHits(overlaps[i])
    subject <- subjectHits(overlaps[i])
    
    if (query != subject) {
      q_start <- start(gr[query])
      q_end <- end(gr[query])
      s_start <- start(gr[subject])
      s_end <- end(gr[subject])
      
      if (q_start < s_start && q_end > s_start) {
        end(gr[query]) <- s_start - 1
      }
      
      if (s_start < q_start && s_end > q_start) {
        end(gr[subject]) <- q_start - 1
      }
    }
  }
  
  gr <- gr[start(gr) <= end(gr)]
  
  return(gr)
}

write_bed_file <- function(gr, output_file, headers) {
  bed_out <- data.frame(seqnames = seqnames(gr),
                        chromStart = start(gr),
                        chromEnd = end(gr),
                        source = mcols(gr)$source,
                        feature = mcols(gr)$feature,
                        score = mcols(gr)$score,
                        strand = strand(gr),
                        frame = mcols(gr)$frame,
                        attribute = mcols(gr)$attribute)
  
  headers[2:3] <- c("chromStart", "chromEnd")
  
  cat(paste(headers, collapse = "\t"), "\n", file = output_file)
  write.table(bed_out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}


bed_info <- read_bed_file(input_file)
bed <- bed_info$data
headers <- bed_info$headers

gr_trimmed <- trim_overlaps(bed)
write_bed_file(gr_trimmed, output_file, headers)

cat("Trimmed BED file has been written to", output_file, "\n")








# Filter input bed file to include first 3 columns and two values in column 9

# Load necessary library
library(dplyr)
library(stringr)

# Define the file path (adjust as necessary)
file_path <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"

# Read the BED file into a data frame
# Specify only the first 9 columns, assuming there are no headers
bed_data <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", colClasses = "character")

# Function to extract ENST and gene_name values from column 9
extract_values <- function(col9) {
  enst_value <- str_extract(col9, "ENST\\S+")
  gene_name_value <- str_extract(col9, "(?<=gene_name )\\S+")
  # Remove trailing semicolon from both values if they exist
  enst_value <- str_replace(enst_value, ";$", "")
  gene_name_value <- str_replace(gene_name_value, ";$", "")
  return(c(enst_value, gene_name_value))
}

# Apply the function to column 9 and create new columns
values_extracted <- t(apply(bed_data, 1, function(row) extract_values(row[9])))
bed_data <- cbind(bed_data, values_extracted)
colnames(bed_data)[(ncol(bed_data)-1):ncol(bed_data)] <- c("ENST", "gene_name")

# Filter the data to include only the first 3 columns and the extracted ENST and gene_name columns
filtered_data <- bed_data %>%
  select(V1, V2, V3, ENST, gene_name) %>%
  filter(!is.na(ENST) & !is.na(gene_name))

# Write the filtered data to a new BED file
write.table(filtered_data, "filtered_Homo_sapiens.GRCh37.87.bed", 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Print a message indicating the script has finished
cat("Filtering complete. The filtered file has been saved as 'filtered_Homo_sapiens.GRCh37.87.bed'.\n")







# Correct input files

# filtered_Homo_sapiens.GRCh37.87.bed

# sorted_homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_pathogenic.bed

# Sort the region bed file
sort -k1,1 -k2,2n filtered_Homo_sapiens.GRCh37.87.bed > sorted_Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_bedtools_non_scientific.bed





input_file <- "sorted_Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions_bedtools_non_scientific.bed"
output_file <- "filtered_sorted_bed_file.bed"

# Remove all rows with first column value that isn't numeric and sort the file

# Load necessary library
library(dplyr)

# Function to check if a value is numeric or "MT", "X", "Y" without coercion
is_strictly_numeric_or_chr <- function(x) {
  grepl("^\\d+$", x) || x %in% c("MT", "X", "Y")
}

# Read the BED file
bed_file <- input_file
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)

# Filter rows where the first column is strictly numeric or "MT", "X", "Y"
filtered_data <- bed_data %>%
  filter(sapply(V1, is_strictly_numeric_or_chr))

# Convert first two columns to appropriate types for sorting
filtered_data <- filtered_data %>%
  mutate(V1 = factor(V1, levels = c(as.character(1:22), "X", "Y", "MT")), V2 = as.numeric(V2))

# Sort the data by the first two columns
sorted_data <- filtered_data %>%
  arrange(V1, V2)

# Write the sorted data back to a new BED file
write.table(sorted_data, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print("Filtering and sorting complete. The new file is saved as 'filtered_sorted_bed_file.bed'.")







# Remove scientific notation

# Load necessary library
library(dplyr)
library(readr)

# Define the file path (adjust as necessary)
file_path <- "filtered_sorted_bed_file.bed"
output_path <- "filtered_no_scientific_notation.bed"

# Read the BED file into a data frame
bed_data <- read_delim(file_path, delim = "\t", col_names = FALSE)

# Function to check for scientific notation in any column
contains_scientific_notation <- function(row) {
  any(grepl("e\\+", row))
}

# Filter out rows that contain scientific notation
filtered_data <- bed_data %>% 
  filter(!apply(., 1, contains_scientific_notation))

# Write the filtered data to a new BED file
write_delim(filtered_data, output_path, delim = "\t", col_names = FALSE)

# Print a message indicating the script has finished
cat("Filtering complete. The filtered file has been saved as 'filtered_no_scientific_notation.bed'.\n")
