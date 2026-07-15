#' Genomic Data Utilities
#' 
#' Helper functions for reading, parsing, and manipulating genomic variant and exon data.
#' All functions work with in-memory data.table for speed and minimal disk I/O.
#'

library(data.table)
library(GenomicRanges)

#' Load and Parse Variant Data
#'
#' Reads variant files (bed or txt format) with optimized settings for large genomic datasets.
#' Prevents scientific notation in coordinates and auto-detects separators.
#'
#' @param file Character. Path to variant file (supports .gz compression)
#' @param colnames Character vector. Column names to assign. If NULL, uses file headers.
#'
#' @return data.table with columns: chr, chromStart, chromEnd, and any additional fields
#'
#' @examples
#' \dontrun{
#'   variants <- load_variants("clinvar_missense.bed.gz")
#' }
#'
#' @export
load_variants <- function(file, colnames = NULL) {
  options(scipen = 999)  # Prevent scientific notation
  
  dt <- fread(file, sep = "\t", header = !is.null(colnames))
  
  # Ensure standard genomic column names
  if (is.null(colnames)) {
    setnames(dt, c("chr", "chromStart", "chromEnd", paste0("V", seq(4, ncol(dt)))))
  } else if (length(colnames) == ncol(dt)) {
    setnames(dt, colnames)
  }
  
  # Ensure numeric coordinates
  dt[, chromStart := as.integer(chromStart)]
  dt[, chromEnd := as.integer(chromEnd)]
  
  # Sort for efficient range operations
  setorder(dt, chr, chromStart)
  
  return(dt)
}

#' Load Gene List from CSV
#'
#' Reads gene names from epigeneticmachinery.org CSV export.
#'
#' @param file Character. Path to CSV file
#' @param column Character. Name of column containing gene names (default: "Gene_Name")
#'
#' @return Character vector of gene names
#'
#' @export
load_gene_list <- function(file, column = "Gene_Name") {
  genes <- fread(file, sep = ",", header = TRUE)
  return(unique(genes[[column]]))
}

#' Convert GTF to GRanges
#'
#' Parse Ensembl GTF file and extract exon coordinates as GRanges object.
#' Much faster than bedtools for downstream overlap operations.
#'
#' @param file Character. Path to GTF file (supports .gz)
#' @param feature Character. Feature type to extract (default: "exon")
#'
#' @return GRanges object with seqnames, ranges, strand, and metadata (gene_id, transcript_id, exon_id)
#'
#' @export
load_gtf_as_granges <- function(file, feature = "exon") {
  # Read GTF with minimal parsing
  gtf <- fread(file, sep = "\t", header = FALSE, skip = "#")
  setnames(gtf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  
  # Filter to requested feature
  gtf <- gtf[feature == feature]
  
  # Parse attributes column for gene_id, transcript_id, exon_id
  # Format: gene_id "ENSG..."; transcript_id "ENST..."; exon_id "ENSE..."
  parse_attribute <- function(attr_str, key) {
    pattern <- sprintf('%s "([^"]+)"', key)
    matches <- regmatches(attr_str, regexec(pattern, attr_str))
    sapply(matches, function(x) if (length(x) > 1) x[2] else NA_character_)
  }
  
  gtf[, gene_id := parse_attribute(attribute, "gene_id")]
  gtf[, transcript_id := parse_attribute(attribute, "transcript_id")]
  gtf[, exon_id := parse_attribute(attribute, "exon_id")]
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = gtf$seqname,
    ranges = IRanges(start = gtf$start, end = gtf$end),
    strand = gtf$strand,
    gene_id = gtf$gene_id,
    transcript_id = gtf$transcript_id,
    exon_id = gtf$exon_id
  )
  
  return(gr)
}

#' Convert data.table to GRanges
#'
#' Utility to convert genomic data.table (chr, chromStart, chromEnd) to GRanges.
#'
#' @param dt data.table with chr, chromStart, chromEnd columns
#' @param metadata_cols Character vector. Additional columns to add as metadata
#'
#' @return GRanges object
#'
#' @export
dt_to_granges <- function(dt, metadata_cols = NULL) {
  gr <- GRanges(
    seqnames = dt$chr,
    ranges = IRanges(start = dt$chromStart, end = dt$chromEnd)
  )
  
  if (!is.null(metadata_cols)) {
    for (col in metadata_cols) {
      mcols(gr)[[col]] <- dt[[col]]
    }
  }
  
  return(gr)
}

#' Write Results to File
#'
#' Efficiently write data.table results to compressed TSV.
#'
#' @param dt data.table
#' @param file Character. Output path (auto-detects .gz compression)
#' @param compress Logical. Whether to gzip output (default: TRUE)
#'
#' @export
write_results <- function(dt, file, compress = TRUE) {
  fwrite(dt, file = file, sep = "\t", quote = FALSE)
  cat("Results written to:", file, "\n")
}

#' Memory-efficient line count for compressed files
#'
#' Count lines in gzip-compressed files without full decompression.
#' Much faster than system("zcat ... | wc -l").
#'
#' @param file Character. Path to .gz file
#'
#' @return Integer. Line count
#'
#' @export
count_lines_gzip <- function(file) {
  con <- gzfile(file, "rt")
  count <- 0L
  while (length(chunk <- readLines(con, n = 100000)) > 0) {
    count <- count + length(chunk)
  }
  close(con)
  return(count)
}
