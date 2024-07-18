# Main script
input_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"

# Load necessary libraries
library(dplyr)
library(GenomicRanges)

# Function to read and parse the BED file
read_bed_file <- function(file) {
  bed <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end", "source", "feature", "score", "strand", "frame", "attribute")
  
  # Ensure start and end columns are numeric
  bed$start <- suppressWarnings(as.numeric(bed$start))
  bed$end <- suppressWarnings(as.numeric(bed$end))
  
  # Filter out rows with NA values in start or end after coercion
  bed <- bed %>% filter(!is.na(start) & !is.na(end))
  
  return(list(data = bed, headers = colnames(bed)))
}

# Function to trim overlapping regions
trim_overlaps <- function(bed) {
  gr <- GRanges(seqnames = bed$chr,
                ranges = IRanges(start = bed$start, end = bed$end),
                strand = bed$strand,
                source = bed$source,
                feature = bed$feature,
                score = bed$score,
                frame = bed$frame,
                attribute = bed$attribute)
  
  # Reduce overlapping ranges
  reduced_gr <- reduce(gr, min.gapwidth = 0)
  
  # Trim the original GRanges based on the reduced ranges
  trimmed_gr <- gr
  
  for (i in seq_along(reduced_gr)) {
    overlaps <- findOverlaps(trimmed_gr, reduced_gr[i])
    if (length(overlaps) > 1) {
      overlap_hits <- subjectHits(overlaps)
      start_trim <- start(reduced_gr[i])
      end_trim <- end(reduced_gr[i])
      
      for (hit in overlap_hits) {
        if (start(trimmed_gr[hit]) < start_trim) {
          end(trimmed_gr[hit]) <- start_trim - 1
        }
        if (end(trimmed_gr[hit]) > end_trim) {
          start(trimmed_gr[hit]) <- end_trim + 1
        }
      }
    }
  }
  
  # Ensure no negative lengths
  trimmed_gr <- trimmed_gr[start(trimmed_gr) <= end(trimmed_gr)]
  
  return(trimmed_gr)
}

# Function to write the trimmed BED file with headers
write_bed_file <- function(gr, output_file, headers) {
  bed_out <- data.frame(seqnames = seqnames(gr),
                        start = start(gr),
                        end = end(gr),
                        source = mcols(gr)$source,
                        feature = mcols(gr)$feature,
                        score = mcols(gr)$score,
                        strand = strand(gr),
                        frame = mcols(gr)$frame,
                        attribute = mcols(gr)$attribute)
  # Write headers
  cat(paste(headers, collapse = "\t"), "\n", file = output_file)
  # Write data
  write.table(bed_out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# Read input BED file, limited to the first 1000 lines
bed_info <- read_bed_file(input_file)
bed <- bed_info$data
headers <- bed_info$headers

# Trim overlapping regions
gr_trimmed <- trim_overlaps(bed)

# Write trimmed BED file
write_bed_file(gr_trimmed, output_file, headers)

cat("Variant counts have been added and data have been written to", output_file, "\n")
