# Main script
input_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff.bed"
output_file <- "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff_isoform_specific_regions.bed"

# Load necessary libraries
library(dplyr)
library(GenomicRanges)

# Function to read and parse the BED file
read_bed_file <- function(file) {
  bed <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end", "source", "feature", "score", "strand", "frame", "attribute")
  return(bed)
}

# Function to check for overlaps and trim regions
trim_overlaps <- function(bed) {
  gr <- GRanges(seqnames = bed$chr,
                ranges = IRanges(start = bed$start, end = bed$end),
                strand = bed$strand,
                score = bed$score,
                attribute = bed$attribute)
  
  # Find overlaps
  overlaps <- findOverlaps(gr, gr)
  
  # Initialize a list to store trimmed regions
  trimmed_regions <- list()
  
  # Process each overlap
  for (i in seq_along(overlaps)) {
    overlap <- overlaps[i]
    if (queryHits(overlap) != subjectHits(overlap)) {
      q_start <- start(gr[queryHits(overlap)])
      q_end <- end(gr[queryHits(overlap)])
      s_start <- start(gr[subjectHits(overlap)])
      s_end <- end(gr[subjectHits(overlap)])
      
      if (q_start < s_start & q_end > s_start) {
        end(gr[queryHits(overlap)]) <- s_start - 1
      }
      
      if (s_start < q_start & s_end > q_start) {
        end(gr[subjectHits(overlap)]) <- q_start - 1
      }
    }
  }
  
  return(gr)
}

# Function to write the trimmed BED file
write_bed_file <- function(gr, output_file) {
  bed_out <- data.frame(seqnames = seqnames(gr),
                        start = start(gr),
                        end = end(gr),
                        source = mcols(gr)$source,
                        feature = mcols(gr)$feature,
                        score = mcols(gr)$score,
                        strand = strand(gr),
                        frame = mcols(gr)$frame,
                        attribute = mcols(gr)$attribute)
  write.table(bed_out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Read input BED file
bed <- read_bed_file(input_file)

# Trim overlapping regions
gr_trimmed <- trim_overlaps(bed)

# Write trimmed BED file
write_bed_file(gr_trimmed, output_file)
