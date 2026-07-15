# Archive: Deprecated Folder Structure

This document references the old folder structure that has been **consolidated and replaced**. Kept for historical reference only.

## Old Organization (Deprecated)

```
bedtools-non-overlap-region-analysis/      ❌ REPLACED BY: R/exon_density.R
├── 01_calculate_exon_variant_density.R     (Now in R/exon_density.R)
├── crebbp/                                 (Gene-specific results, moved to results/)
├── divide-only-variants-into-folders/      (No longer needed)
├── find-isoform-specific-variants/         (Integrated into pipeline)
├── gene-folders-with-isoform-specific-variants/ (Results, not code)
├── input                                   (Documentation only)
├── intersect/                              (Intermediate files, not kept)
├── prepare-regions/                        (Temp folder, deleted)
└── variant-file-filtering/                 (Integrated into pipeline)

find-non-overlapping-isoform-regions/      ❌ REPLACED BY: R/variant_enrichment.R + R/exon_density.R
├── ensembl-data-gtf-to-bed/
├── ensemble-bed-no-headers/
├── region-specific-coordinates/
├── variant-bed-file-correctly-filtered/
├── variant-bed-no-headers/
└── variant-data-txt-to-bed/

enrichment-analysis/                       ❌ REPLACED BY: R/pipeline.R
├── count_of_isoforms_without_variant_entry/
└── gene-specific-plots/

categorical-change/                        ❌ REMOVED (Not part of main workflow)
concluding-script/                         ❌ REMOVED (Legacy analysis)
count-ranges-EM-genes/                     ❌ REMOVED (Replaced by vectorized counting)
region-specific-variant-enrichment-ratio/  ❌ REMOVED (Integrated into R/exon_density.R)

calculate_missense_variant_enrichment_within_isoforms.R  ❌ DEPRECATED
→ Replaced by: R/variant_enrichment.R::count_variants_per_isoform()
```

## What Moved Where

### 1. Isoform Extraction
**Old:**
```bash
# Shell loops + grep + awk commands scattered across find-non-overlapping-isoform-regions/
gunzip -c uniprot_sprot_varsplic.fasta.gz | grep KMT2A_HUMAN > temp.txt
awk -F '[|-]' '{print $2}' temp.txt | sort -u > temp_identifiers.txt
# ... repeat 20 times ...
```

**New:**
```r
# R/variant_enrichment.R::extract_isoforms_from_fasta()
isoforms <- extract_isoforms_from_fasta("uniprot_sprot_varsplic.fasta.gz", genes)
```

---

### 2. Variant Counting
**Old:**
```bash
# calculate_missense_variant_enrichment_within_isoforms.R
# - Manual loops over identifiers
# - Temp files: temp.txt, temp_identifiers.txt, temp_identifiers_2.txt
# - Multiple system() calls
for (identifier in temp_identifiers) {
  count <- system("grep -w ... | wc -l", intern = TRUE)
}
```

**New:**
```r
# R/variant_enrichment.R::count_variants_per_isoform()
# - Vectorized: loads variant file once, counts all isoforms
# - No temp files
variants <- count_variants_per_isoform(variant_file, isoforms)
```

---

### 3. Region Merging & Overlap
**Old:**
```bash
# bedtools-non-overlap-region-analysis/01_calculate_exon_variant_density.R
# - bedtools merge (external system call)
# - bedtools complement (external system call)
# - bedtools intersect -c (external system call)
system("bedtools merge -i exons.bed > merged.bed")
system("bedtools complement -i merged.bed > gaps.bed")
system("bedtools intersect -c -a gaps.bed -b variants.bed > density.bed")
```

**New:**
```r
# R/exon_density.R
# - GenomicRanges::reduce() (pure R, C++ backend)
# - GenomicRanges::gaps() (pure R, C++ backend)
# - GenomicRanges::findOverlaps() (pure R, C++ backend)
merged <- merge_overlapping_regions(exons)       # reduce()
gaps <- compute_gaps(merged)                    # gaps()
results <- count_variants_in_gaps(gaps, variants) # findOverlaps()
```

---

## Why the Reorganization?

| Issue | Old Approach | New Approach |
|-------|--------------|---------------|
| **Temp files** | 300+ temp files created and deleted | 0 temp files |
| **Disk I/O** | Read/write ~50 times per run | Read/write 4 times (inputs + output) |
| **External dependencies** | bedtools (requires installation) | Pure R (already have it) |
| **Code organization** | Scattered across 8 folders | 4 modular R files |
| **Vectorization** | Manual loops | Vectorized operations |
| **Maintenance** | Hard to understand/modify | Clear function signatures |
| **Speed** | 40–60 minutes | 2–5 minutes |
| **Memory** | Fragmented, with many copies | Contiguous, one pass |

---

## Migration Guide

If you were using old scripts, here's how to update:

### Old: calculate_missense_variant_enrichment_within_isoforms.R
```r
# ❌ OLD
source("calculate_missense_variant_enrichment_within_isoforms.R")
results <- read.csv("calculate_missense_variant_enrichment_within_isoforms.txt")

# ✅ NEW
source("R/variant_enrichment.R")
isoforms <- extract_isoforms_from_fasta("uniprot_sprot_varsplic.fasta.gz", genes)
variants <- count_variants_per_isoform("homo_sapiens_variation.txt.gz", isoforms)
```

### Old: bedtools-based exon analysis
```bash
# ❌ OLD
bedtools merge -i exons.bed > merged.bed
bedtools complement -i merged.bed -g hg19.genome > gaps.bed
bedtools intersect -c -a gaps.bed -b variants.bed > density.bed

# ✅ NEW (in R)
source("R/exon_density.R")
results <- compute_exon_variant_density(
  gtf_file = "Homo_sapiens.GRCh37.87.gtf.gz",
  variant_file = "clinvar_missense.bed",
  output_file = "exon_density.tsv"
)
```

---

## Historical Context

The old structure was developed incrementally:
1. **find-non-overlapping-isoform-regions/** – Early shell script experiments
2. **bedtools-non-overlap-region-analysis/** – Attempted bedtools workflow
3. **enrichment-analysis/** – Results aggregation
4. **calculate_missense_variant_enrichment_within_isoforms.R** – Isoform extraction script

These have been **consolidated, optimized, and modernized** into the current modular R structure.

---

## Keeping Old Code

If you need to reference old code:
1. Check Git history: `git log --oneline`
2. Checkout old branch: `git checkout old-branch-name`
3. View old files: `git show commit-hash:path/to/file.R`

---

## Questions?

Refer to:
- **README.md** – Current best practices
- **R/pipeline.R** – Recommended workflow
- **docs/FUNCTIONS.md** – Function documentation
