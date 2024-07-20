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

# Main script
input_file <- "Homo_sapiens.GRCh37.87_with_headers.bed"
output_file <- "Homo_sapiens.GRCh37.87_with_headers_isoform_specific_regions.bed"

bed_info <- read_bed_file(input_file)
bed <- bed_info$data
headers <- bed_info$headers

gr_trimmed <- trim_overlaps(bed)
write_bed_file(gr_trimmed, output_file, headers)

cat("Trimmed BED file has been written to", output_file, "\n")
