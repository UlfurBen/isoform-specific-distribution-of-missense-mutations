#' Exon Variant Density Analysis
#'
#' Replaces bedtools workflow entirely with GenomicRanges operations.
#' In-memory computation, no temp files, 10-100x faster.
#'

library(data.table)
library(GenomicRanges)

#' Merge Overlapping Regions (Replaces bedtools merge)
#'
#' Combines overlapping exon regions into unified contiguous blocks.
#' Uses GenomicRanges::reduce() for speed and correctness.
#'
#' @param regions GRanges or data.table. Genomic regions (chr, start, end)
#'
#' @return GRanges object with merged, non-overlapping regions
#'
#' @export
merge_overlapping_regions <- function(regions) {
  if (is(regions, "data.table")) {
    regions <- dt_to_granges(regions)
  }
  
  # Reduce overlapping intervals by chromosome and strand
  merged <- reduce(split(regions, seqnames(regions)))
  merged <- unlist(merged)
  
  return(merged)
}

#' Compute Non-Overlapping Gaps (Replaces bedtools complement)
#'
#' Identifies gaps between merged exon regions.
#' These are the "unique exon" regions of interest for mutation enrichment.
#'
#' @param merged_regions GRanges. Output from merge_overlapping_regions()
#' @param genome_lengths Named integer vector. Chromosome lengths (e.g., from seqinfo)
#'
#' @return GRanges object representing gaps between merged regions
#'
#' @export
compute_gaps <- function(merged_regions, genome_lengths = NULL) {
  gaps <- gaps(merged_regions)
  
  # Filter out regions beyond chromosome ends if genome_lengths provided
  if (!is.null(genome_lengths)) {
    gaps <- gaps[end(gaps) <= genome_lengths[as.character(seqnames(gaps))]]
  }
  
  return(gaps)
}

#' Count Variants in Gaps (Replaces bedtools intersect -c)
#'
#' Fast overlap join: counts how many variants fall within each gap region.
#' Uses GenomicRanges::findOverlaps() for accuracy and speed.
#'
#' @param gaps GRanges. Gap regions (output from compute_gaps())
#' @param variants GRanges or data.table. Variant locations
#'
#' @return data.table with columns: chr, start, end, variant_count
#'
#' @export
count_variants_in_gaps <- function(gaps, variants) {
  if (is(variants, "data.table")) {
    variants <- dt_to_granges(variants)
  }
  
  cat("Counting variants in", length(gaps), "gap regions...\n")
  
  # Find which variants overlap each gap
  hits <- findOverlaps(gaps, variants, type = "within")
  
  # Count by gap
  counts <- countOverlaps(gaps, variants, type = "within")
  
  # Convert back to data.table
  result <- data.table(
    chr = as.character(seqnames(gaps)),
    chromStart = start(gaps),
    chromEnd = end(gaps),
    variant_count = counts
  )
  
  result[, region_length := chromEnd - chromStart]
  result[, variant_density := ifelse(region_length > 0, variant_count / region_length, 0)]
  
  return(result[order(-variant_density)])
}

#' Complete Pipeline: GTF + Variants → Density Analysis
#'
#' End-to-end workflow combining region processing and variant counting.
#' No bedtools, no temp files, all in-memory.
#'
#' @param gtf_file Character. Path to Ensembl GTF.gz
#' @param variant_file Character. Path to variant BED file
#' @param output_file Character. Output TSV path
#' @param feature Character. Feature type to extract (default: "exon")
#'
#' @return data.table (invisibly) and writes results to output_file
#'
#' @export
compute_exon_variant_density <- function(gtf_file, variant_file, output_file, feature = "exon") {
  cat("Starting exon variant density computation...\n\n")
  
  # Load data
  cat("1. Loading GTF regions...\n")
  regions <- load_gtf_as_granges(gtf_file, feature = feature)
  cat("   Loaded", length(regions), "regions\n\n")
  
  cat("2. Loading variants...\n")
  variants <- load_variants(variant_file)
  cat("   Loaded", nrow(variants), "variants\n\n")
  
  # Process regions
  cat("3. Merging overlapping exons...\n")
  merged <- merge_overlapping_regions(regions)
  cat("   Merged to", length(merged), "blocks\n\n")
  
  cat("4. Computing gaps (unique exon regions)...\n")
  gaps <- compute_gaps(merged)
  cat("   Found", length(gaps), "gap regions\n\n")
  
  # Count variants
  cat("5. Counting variants in gap regions...\n")
  result <- count_variants_in_gaps(gaps, variants)
  
  # Save
  cat("\n6. Writing results...\n")
  write_results(result, output_file)
  
  cat("\nDensity computation complete!\n")
  cat("Top 10 highest-density regions:\n")
  print(head(result, 10))
  
  invisible(result)
}
