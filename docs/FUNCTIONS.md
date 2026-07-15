# Function Reference

Complete documentation of all R functions in the pipeline.

## R/utils.R

### load_variants(file, colnames = NULL)
Load and parse genomic variant data.

**Parameters:**
- `file` (character): Path to variant file (BED, TXT, supports .gz)
- `colnames` (character vector, optional): Column names to assign

**Returns:** data.table with chr, chromStart, chromEnd, and metadata

**Example:**
```r
variants <- load_variants("variants.bed.gz")
head(variants)
```

---

### load_gene_list(file, column = "Gene_Name")
Load gene names from epigeneticmachinery.org CSV.

**Parameters:**
- `file` (character): Path to CSV file
- `column` (character): Column name containing gene names

**Returns:** Character vector of unique gene names

**Example:**
```r
genes <- load_gene_list("The-Epigenetic-Machinery.csv")
print(length(genes))  # 300 genes
```

---

### load_gtf_as_granges(file, feature = "exon")
Parse Ensembl GTF and extract regions as GRanges.

**Parameters:**
- `file` (character): Path to GTF file (supports .gz)
- `feature` (character): Feature type to extract ("exon", "transcript", "gene")

**Returns:** GRanges object with gene_id, transcript_id, exon_id metadata

**Example:**
```r
exons <- load_gtf_as_granges("Homo_sapiens.GRCh37.87.gtf.gz")
print(length(exons))  # ~1.2 million exons
```

---

### dt_to_granges(dt, metadata_cols = NULL)
Convert data.table to GRanges object.

**Parameters:**
- `dt` (data.table): Must have chr, chromStart, chromEnd columns
- `metadata_cols` (character vector, optional): Additional columns to add as metadata

**Returns:** GRanges object

**Example:**
```r
variants_gr <- dt_to_granges(variants, metadata_cols = c("variant_id", "effect"))
```

---

### write_results(dt, file, compress = TRUE)
Write results to TSV file.

**Parameters:**
- `dt` (data.table): Data to write
- `file` (character): Output path
- `compress` (logical): Auto-detect .gz compression (default: TRUE)

**Example:**
```r
write_results(results, "output.tsv")
```

---

### count_lines_gzip(file)
Fast line count for gzip-compressed files.

**Parameters:**
- `file` (character): Path to .gz file

**Returns:** Integer line count

**Example:**
```r
line_count <- count_lines_gzip("homo_sapiens_variation.txt.gz")
print(line_count)  # ~53 million lines
```

---

## R/variant_enrichment.R

### extract_isoforms_from_fasta(fasta_file, gene_names)
Extract isoform IDs from UniProt FASTA headers.

**Parameters:**
- `fasta_file` (character): Path to UniProt FASTA file (supports .gz)
- `gene_names` (character vector): Gene names to search (with _HUMAN suffix)

**Returns:** data.table with gene_name, isoform_id, is_canonical columns

**Example:**
```r
isoforms <- extract_isoforms_from_fasta(
  "uniprot_sprot_varsplic.fasta.gz",
  c("KMT2A", "AIRE")
)
head(isoforms)
```

---

### count_variants_per_isoform(variant_file, isoforms, variant_type = "missense_variant")
Vectorized counting of variants per isoform.

**Parameters:**
- `variant_file` (character): Path to variation file (supports .gz)
- `isoforms` (data.table): Output from extract_isoforms_from_fasta()
- `variant_type` (character): Type of variant to count

**Returns:** data.table with isoform_id, variant_count, gene_name (sorted by count)

**Example:**
```r
variants <- count_variants_per_isoform(
  "homo_sapiens_variation.txt.gz",
  isoforms
)
head(variants, 10)
```

---

### add_isoform_lengths(isoform_counts, fasta_file)
Merge sequence lengths and compute variant density.

**Parameters:**
- `isoform_counts` (data.table): Output from count_variants_per_isoform()
- `fasta_file` (character): Path to UniProt FASTA (supports .gz)

**Returns:** data.table with added sequence_length and density columns

**Example:**
```r
density <- add_isoform_lengths(isoform_counts, "uniprot_sprot_varsplic.fasta.gz")
head(density)
```

---

### filter_clinvar_missense(input_file, output_file)
Single-pass extraction of ClinVar missense variants.

**Parameters:**
- `input_file` (character): Path to homo_sapiens_variation.txt.gz
- `output_file` (character): Output file path

**Returns:** Integer count of variants written

**Example:**
```r
count <- filter_clinvar_missense(
  "homo_sapiens_variation.txt.gz",
  "homo_sapiens_variation_missense_ClinVar.bed"
)
cat("Extracted", count, "variants\n")
```

---

## R/exon_density.R

### merge_overlapping_regions(regions)
Merge overlapping exon regions (replaces bedtools merge).

**Parameters:**
- `regions` (GRanges or data.table): Genomic regions

**Returns:** GRanges with merged, non-overlapping regions

**Example:**
```r
exons <- load_gtf_as_granges("Homo_sapiens.GRCh37.87.gtf.gz")
merged <- merge_overlapping_regions(exons)
cat("Merged", length(exons), "exons to", length(merged), "blocks\n")
```

---

### compute_gaps(merged_regions, genome_lengths = NULL)
Identify gaps between merged regions (replaces bedtools complement).

**Parameters:**
- `merged_regions` (GRanges): Output from merge_overlapping_regions()
- `genome_lengths` (named integer vector, optional): Chromosome lengths

**Returns:** GRanges of gap regions

**Example:**
```r
gaps <- compute_gaps(merged_regions)
cat("Found", length(gaps), "gap regions\n")
```

---

### count_variants_in_gaps(gaps, variants)
Count variants falling within gap regions (replaces bedtools intersect -c).

**Parameters:**
- `gaps` (GRanges): Gap regions from compute_gaps()
- `variants` (GRanges or data.table): Variant locations

**Returns:** data.table with chr, chromStart, chromEnd, variant_count, region_length, variant_density

**Example:**
```r
results <- count_variants_in_gaps(gaps, variants_gr)
head(results)
```

---

### compute_exon_variant_density(gtf_file, variant_file, output_file, feature = "exon")
Complete end-to-end exon variant density pipeline.

**Parameters:**
- `gtf_file` (character): Path to Ensembl GTF.gz
- `variant_file` (character): Path to variant BED file
- `output_file` (character): Output TSV path
- `feature` (character): Feature type to extract

**Returns:** data.table (invisibly), writes results to output_file

**Example:**
```r
results <- compute_exon_variant_density(
  "Homo_sapiens.GRCh37.87.gtf.gz",
  "clinvar_missense.bed",
  "exon_density.tsv"
)
head(results, 10)
```

---

## R/pipeline.R

### run_isoform_enrichment_pipeline(...)
Complete workflow orchestrator (main entry point).

**Parameters:**
- `gene_list_file` (character): Path to The-Epigenetic-Machinery.csv
- `uniprot_fasta` (character): Path to uniprot_sprot_varsplic.fasta.gz
- `uniprot_variants` (character): Path to homo_sapiens_variation.txt.gz
- `gtf_file` (character): Path to Homo_sapiens.GRCh37.87.gtf.gz
- `clinvar_variants` (character): Path to ClinVar missense BED
- `output_dir` (character): Output directory (default: "./results")
- `n_jobs` (integer): Parallel chunks for processing (default: 4)

**Returns:** Invisible list with all intermediate results

**Example:**
```r
results <- run_isoform_enrichment_pipeline(
  gene_list_file = "The-Epigenetic-Machinery.csv",
  uniprot_fasta = "uniprot_sprot_varsplic.fasta.gz",
  uniprot_variants = "homo_sapiens_variation.txt.gz",
  gtf_file = "Homo_sapiens.GRCh37.87.gtf.gz",
  clinvar_variants = "clinvar_missense.bed",
  output_dir = "./results"
)

# Access individual results
print(results$isoforms)
print(results$exon_density)
```

---

### quick_start()
Run pipeline with pre-configured paths for standard data locations.

**Parameters:** None

**Returns:** Same as run_isoform_enrichment_pipeline()

**Example:**
```r
quick_start()  # Assumes data files in working directory
```

---

## Performance Notes

### Vectorization
All functions use vectorized operations. Avoid manual loops:

```r
# ❌ SLOW: Manual loop
for (i in 1:nrow(variants)) {
  count <- grep(variants$id[i], variation_db)
}

# ✅ FAST: Vectorized
count_variants_per_isoform(variation_db, variants)  # Loads entire file, counts all at once
```

### Memory Efficiency
- Data.table operations are "in-place" (don't copy entire dataset)
- GRanges overlap operations use efficient C++ backends
- Streaming for very large files (future enhancement)

---
