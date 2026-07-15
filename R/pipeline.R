#' Main Analysis Pipeline
#'
#' Orchestrates the complete workflow:
#' 1. Load gene list from epigeneticmachinery.org CSV
#' 2. Extract isoforms from UniProt FASTA
#' 3. Count missense variants per isoform (ClinVar + all DBs)
#' 4. Compute exon variant density (replaces bedtools)
#' 5. Generate summary statistics
#'
#' All operations are in-memory and vectorized for speed.
#' No system() calls, no temp files, no bedtools dependency.
#'

library(data.table)
library(GenomicRanges)

#' Run Complete Variant Enrichment Pipeline
#'
#' @param gene_list_file Character. Path to The-Epigenetic-Machinery.csv
#' @param uniprot_fasta Character. Path to uniprot_sprot_varsplic.fasta.gz
#' @param uniprot_variants Character. Path to homo_sapiens_variation.txt.gz
#' @param gtf_file Character. Path to Homo_sapiens.GRCh37.87.gtf.gz
#' @param clinvar_variants Character. Path to clinvar missense BED/txt file
#' @param output_dir Character. Directory for all outputs
#' @param n_jobs Integer. Parallel chunks for variant processing (default: 4)
#'
#' @return Invisible list with all intermediate results and file paths
#'
#' @export
run_isoform_enrichment_pipeline <- function(
    gene_list_file,
    uniprot_fasta,
    uniprot_variants,
    gtf_file,
    clinvar_variants,
    output_dir = "./results",
    n_jobs = 4) {
  
  start_time <- Sys.time()
  
  # Setup
  dir.create(output_dir, showWarnings = FALSE)
  cat("========================================\n")
  cat("Isoform Variant Enrichment Pipeline\n")
  cat("========================================\n\n")
  
  # Step 1: Load genes
  cat("STEP 1: Loading gene list...\n")
  genes <- load_gene_list(gene_list_file)
  cat("Loaded", length(genes), "genes from", basename(gene_list_file), "\n\n")
  
  # Step 2: Extract isoforms
  cat("STEP 2: Extracting isoforms from UniProt FASTA...\n")
  isoforms <- extract_isoforms_from_fasta(uniprot_fasta, genes)
  cat("Found", nrow(isoforms), "isoforms across", n_distinct(isoforms$gene_name), "genes\n")
  cat("Canonical:", sum(isoforms$is_canonical), "|  Splice variants:", sum(!isoforms$is_canonical), "\n\n")
  
  # Step 3: Count variants per isoform
  cat("STEP 3: Counting missense variants per isoform...\n")
  isoform_counts <- count_variants_per_isoform(uniprot_variants, isoforms)
  cat("Mean variants per isoform:", round(mean(isoform_counts$variant_count), 1), "\n")
  cat("Median variants per isoform:", median(isoform_counts$variant_count), "\n")
  cat("Max variants per isoform:", max(isoform_counts$variant_count), "\n\n")
  
  # Step 4: Add sequence lengths and compute density
  cat("STEP 4: Computing variant density (variants per bp)...\n")
  isoform_density <- add_isoform_lengths(isoform_counts, uniprot_fasta)
  cat("Mean density:", round(mean(isoform_density$density, na.rm = TRUE), 4), "variants/bp\n")
  cat("Isoforms with density > 0.01:", sum(isoform_density$density > 0.01, na.rm = TRUE), "\n\n")
  
  # Step 5: Exon-level analysis (replaces bedtools)
  cat("STEP 5: Computing exon variant density (GenomicRanges, no bedtools)...\n")
  exon_density <- compute_exon_variant_density(
    gtf_file = gtf_file,
    variant_file = clinvar_variants,
    output_file = file.path(output_dir, "exon_variant_density.tsv")
  )
  cat("\n")
  
  # Save intermediate results
  cat("STEP 6: Writing results to disk...\n")
  write_results(isoforms, file.path(output_dir, "01_isoforms.tsv"))
  write_results(isoform_counts, file.path(output_dir, "02_variant_counts.tsv"))
  write_results(isoform_density, file.path(output_dir, "03_variant_density.tsv"))
  
  # Summary statistics
  cat("\nSTEP 7: Summary Statistics\n")
  cat("============================\n")
  cat("Total genes:", n_distinct(isoforms$gene_name), "\n")
  cat("Total isoforms:", nrow(isoforms), "\n")
  cat("Total variants:", sum(isoform_counts$variant_count), "\n")
  cat("Total gap regions:", nrow(exon_density), "\n")
  cat("\nHighest-density isoforms:\n")
  print(head(isoform_density[, .(gene_name, isoform_id, variant_count, sequence_length, density)], 10))
  cat("\nHighest-density exon gaps:\n")
  print(head(exon_density, 10))
  
  elapsed <- Sys.time() - start_time
  cat("\n========================================\n")
  cat("Pipeline completed in", round(elapsed[[1]], 1), "seconds\n")
  cat("Results saved to:", output_dir, "\n")
  cat("========================================\n")
  
  return(invisible(list(
    genes = genes,
    isoforms = isoforms,
    isoform_counts = isoform_counts,
    isoform_density = isoform_density,
    exon_density = exon_density,
    output_dir = output_dir
  )))
}

#' Quick Start: Run with Default Paths
#'
#' Assumes standard data file locations in current directory.
#' Edit paths as needed for your setup.
#'
#' @export
quick_start <- function() {
  run_isoform_enrichment_pipeline(
    gene_list_file = "The-Epigenetic-Machinery.csv",
    uniprot_fasta = "uniprot_sprot_varsplic.fasta.gz",
    uniprot_variants = "homo_sapiens_variation.txt.gz",
    gtf_file = "Homo_sapiens.GRCh37.87.gtf.gz",
    clinvar_variants = "homo_sapiens_variation_missense_ClinVar.bed",
    output_dir = "./results"
  )
}
