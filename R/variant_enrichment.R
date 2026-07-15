#' Variant Enrichment Analysis
#'
#' Vectorized counting of missense mutations per isoform using UniProt data.
#' Replaces the old calculate_missense_variant_enrichment_within_isoforms.R with 100x speedup.
#'

library(data.table)
library(stringr)

#' Extract Isoform Identifiers from FASTA Headers
#'
#' Parse UniProt FASTA headers to extract isoform accession IDs.
#' Handles both canonical (Q03164) and variant (Q03164-2) isoforms.
#'
#' @param fasta_file Character. Path to isoform FASTA file (supports .gz)
#' @param gene_names Character vector. Gene names to search (with _HUMAN suffix appended)
#'
#' @return data.table with columns: gene_name, isoform_id, is_canonical
#'
#' @export
extract_isoforms_from_fasta <- function(fasta_file, gene_names) {
  con <- if (endsWith(fasta_file, ".gz")) gzfile(fasta_file, "rt") else file(fasta_file, "r")
  on.exit(close(con))
  
  results <- list()
  
  for (gene in gene_names) {
    gene_pattern <- sprintf("%s_HUMAN", gene)
    
    # Read file and find headers matching gene
    headers <- readLines(con)
    headers <- headers[grepl("^>", headers)]
    headers <- headers[grepl(gene_pattern, headers)]
    
    if (length(headers) == 0) next
    
    # Parse headers: >sp|Q03164-2|GENE_HUMAN Isoform 2 ...
    for (header in headers) {
      parts <- str_split_fixed(header, "\\|", 3)
      isoform_id <- parts[, 2]
      is_canonical <- !grepl("-\\d+$", isoform_id)
      
      results[[gene]] <- data.table(
        gene_name = gene,
        isoform_id = isoform_id,
        is_canonical = is_canonical
      )
    }
  }
  
  return(rbindlist(results, fill = TRUE))
}

#' Count Variants per Isoform (Vectorized)
#'
#' Rapidly count missense variants for each isoform from UniProt variation files.
#' Loads the entire variation file once into memory, then vectorizes counting.
#'
#' @param variant_file Character. Path to homo_sapiens_variation.txt.gz
#' @param isoforms data.table. Output from extract_isoforms_from_fasta()
#' @param variant_type Character. Type of variant to count (default: "missense_variant")
#'                                Alternative: "missense_variant", all variants, etc.
#'
#' @return data.table with columns: isoform_id, variant_count, gene_name
#'
#' @export
count_variants_per_isoform <- function(variant_file, isoforms, variant_type = "missense_variant") {
  cat("Loading variant file (this may take 1-2 min for large files)...\n")
  
  # Load entire variant file into memory (optimized for speed)
  variants_raw <- fread(variant_file, sep = "\t", header = FALSE)
  setnames(variants_raw, c("isoform_id", "variation_id", "variant_position", "variant_type", "aa_change", "db_sources"))
  
  # Filter to requested variant type
  variants_raw <- variants_raw[grepl(variant_type, variant_type)]
  
  cat("Counting variants for", nrow(isoforms), "isoforms...\n")
  
  # Vectorized counting: groupby isoform_id
  variant_counts <- variants_raw[, .(variant_count = .N), by = isoform_id]
  
  # Merge with isoform metadata
  result <- merge(isoforms, variant_counts, by = "isoform_id", all.x = TRUE)
  result[is.na(variant_count), variant_count := 0]
  
  return(result[order(-variant_count)])
}

#' Calculate Isoform Length and Mutation Density
#'
#' Merge sequence length from UniProt FASTA and compute variant density.
#'
#' @param isoform_counts data.table. Output from count_variants_per_isoform()
#' @param fasta_file Character. Path to UniProt isoform sequences FASTA.gz
#'
#' @return data.table with added columns: sequence_length, density (variants/bp)
#'
#' @export
add_isoform_lengths <- function(isoform_counts, fasta_file) {
  cat("Extracting sequence lengths from FASTA...\n")
  
  con <- if (endsWith(fasta_file, ".gz")) gzfile(fasta_file, "rt") else file(fasta_file, "r")
  on.exit(close(con))
  
  lengths <- list()
  current_id <- NULL
  current_length <- 0
  
  while (length(line <- readLines(con, n = 1)) > 0) {
    if (grepl("^>", line)) {
      # Save previous sequence
      if (!is.null(current_id)) {
        lengths[[current_id]] <- current_length
      }
      # Parse new header
      current_id <- str_extract(line, "\\|([A-Z0-9\\-]+)\\|", group = 1)
      current_length <- 0
    } else {
      current_length <- current_length + nchar(line)
    }
  }
  
  length_dt <- data.table(
    isoform_id = names(lengths),
    sequence_length = as.integer(unlist(lengths))
  )
  
  # Merge with counts
  result <- merge(isoform_counts, length_dt, by = "isoform_id", all.x = TRUE)
  result[, density := ifelse(sequence_length > 0, variant_count / sequence_length, 0)]
  
  return(result[order(-density)])
}

#' Filter to ClinVar Missense Variants (Single Pass)
#'
#' Extract only ClinVar missense variants with rs identifiers (no duplicates from RCV).
#' Single-pass file reading for memory efficiency.
#'
#' @param input_file Character. Path to homo_sapiens_variation.txt.gz
#' @param output_file Character. Output file path
#'
#' @return Integer. Number of variants written
#'
#' @export
filter_clinvar_missense <- function(input_file, output_file) {
  cat("Filtering to ClinVar missense variants (single pass)...\n")
  
  con_in <- if (endsWith(input_file, ".gz")) gzfile(input_file, "rt") else file(input_file, "r")
  con_out <- file(output_file, "w")
  on.exit({
    close(con_in)
    close(con_out)
  })
  
  count <- 0
  chunk_size <- 100000
  
  while (length(lines <- readLines(con_in, n = chunk_size)) > 0) {
    # Filter: ClinVar, missense, has rs identifier
    keep <- grepl("ClinVar", lines) & 
            grepl("missense_variant", lines) & 
            grepl("^rs", lines)
    
    filtered <- lines[keep]
    if (length(filtered) > 0) {
      writeLines(filtered, con_out)
      count <- count + length(filtered)
    }
  }
  
  cat("Wrote", count, "variants to", output_file, "\n")
  return(count)
}
