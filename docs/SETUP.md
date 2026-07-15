# Setup Guide

## System Requirements

- **R** ≥ 4.0 (tested on 4.2+)
- **RAM** ≥ 8 GB (for loading full variant files)
- **Disk** ≥ 10 GB (for data files + results)
- **OS** Linux, macOS, or Windows (with WSL recommended for Windows)

## Step-by-Step Installation

### 1. Clone Repository

```bash
git clone https://github.com/UlfurBen/isoform-specific-distribution-of-missense-mutations.git
cd isoform-specific-distribution-of-missense-mutations
```

### 2. Install R Dependencies

Launch R and run:

```r
# Install BiocManager (if not already installed)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "GenomicRanges",  # Core genomic interval operations
    "IRanges"         # Fast interval arithmetic
))

# Install CRAN packages
install.packages(c(
    "data.table",     # Fast data frame operations
    "stringr"         # String manipulation
))

# Verify installation
library(GenomicRanges)
library(data.table)
library(stringr)
cat("All packages loaded successfully!\n")
```

### 3. Download Data Files

Create a `data/` directory and download required files:

```bash
mkdir -p data
cd data

# UniProt FASTA (1.5 GB) - takes ~10 minutes
echo "Downloading UniProt isoforms..."
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/uniprot_sprot_varsplic.fasta.gz &

# UniProt Variants (4 GB) - takes ~20 minutes
echo "Downloading UniProt variants..."
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/homo_sapiens_variation.txt.gz &

# Ensembl GTF (800 MB)
echo "Downloading Ensembl GTF..."
wget https://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz &

wait  # Wait for all downloads to complete
cd ..
```

**Alternative:** Use `curl` instead of `wget`:
```bash
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/uniprot_sprot_varsplic.fasta.gz
```

### 4. Get Gene List (Manual)

1. Visit https://www.epigeneticmachinery.org/
2. Download "The Epigenetic Machinery" CSV file
3. Save as `The-Epigenetic-Machinery.csv` in repo root

Or programmatically (if available):
```bash
cd data
wget https://www.epigeneticmachinery.org/download/gene_list.csv -O The-Epigenetic-Machinery.csv
cd ..
```

### 5. Create ClinVar Missense Variant File

Optional: Extract only ClinVar missense variants for faster processing:

```bash
# Filter to ClinVar missense variants
zcat data/homo_sapiens_variation.txt.gz | \
  grep -E "^rs" | \
  grep "ClinVar" | \
  grep "missense_variant" > \
  data/homo_sapiens_variation_missense_ClinVar.bed
```

Or in R:
```r
source("R/variant_enrichment.R")
filter_clinvar_missense(
  input_file = "data/homo_sapiens_variation.txt.gz",
  output_file = "data/homo_sapiens_variation_missense_ClinVar.bed"
)
```

### 6. Verify Installation

Run a quick test:

```r
setwd("/path/to/repo")
source("R/utils.R")

# Test loading a small file
cat("Testing file I/O...\n")
genes <- fread("The-Epigenetic-Machinery.csv")
cat("✓ Loaded", nrow(genes), "genes\n")

cat("Installation complete!\n")
```

---

## Configuration

### Custom Data Paths

If your data files are in different locations, modify the paths in `pipeline.R`:

```r
run_isoform_enrichment_pipeline(
  gene_list_file = "/path/to/genes.csv",
  uniprot_fasta = "/path/to/isoforms.fasta.gz",
  uniprot_variants = "/path/to/variants.txt.gz",
  gtf_file = "/path/to/annotations.gtf.gz",
  clinvar_variants = "/path/to/clinvar.bed",
  output_dir = "/path/to/results"
)
```

### Memory Optimization

For systems with <8 GB RAM, pre-filter variants:

```bash
# Keep only missense variants for faster loading
zcat homo_sapiens_variation.txt.gz | \
  grep "missense_variant" | \
  gzip > homo_sapiens_variation_missense_only.txt.gz
```

Then use the smaller file:
```r
run_isoform_enrichment_pipeline(
  ...,
  uniprot_variants = "homo_sapiens_variation_missense_only.txt.gz",
  ...
)
```

---

## Troubleshooting

### "Package not found" Error

```r
# Check if package is installed
require("GenomicRanges")  # FALSE = not installed

# Install from Bioconductor
BiocManager::install("GenomicRanges")
```

### "file not found" Error

```bash
# Check file locations
ls -lh The-Epigenetic-Machinery.csv
ls -lh data/*.gz

# Use absolute paths
getwd()  # See current directory
```

### Out of Memory Error

1. **Reduce variant file size:**
   ```bash
   zcat homo_sapiens_variation.txt.gz | head -1000000 | gzip > variants_small.txt.gz
   ```

2. **Increase R memory limit** (use with caution):
   ```r
   memory.limit(size = 16000)  # Windows only
   ```

3. **Use streaming** (future enhancement)

### Slow Performance

Check RAM availability:
```bash
free -h          # Linux
mem               # macOS
```

Monitor during execution:
```bash
top -p $(pgrep -f Rscript)  # Linux
```

---

## Next Steps

See [README.md](../README.md) for:
- Quick start guide
- Running the pipeline
- Understanding output files
- Common tasks
