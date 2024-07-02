# Define the file paths
gene_csv_file <- "The-Epigenetic-Machinery.csv"
variation_file <- "homo_sapiens_variation.txt.gz"
output_file <- "homo_sapiens_variation_missense_EM_genes.txt"
phenotype_output_file <- "homo_sapiens_variation_missense_EM_genes_phenotype.txt"

# Read gene names from the CSV file
gene_data <- read.csv(gene_csv_file, header = TRUE, stringsAsFactors = FALSE)
gene_names <- gene_data$Gene_Name

# Read the lines of the variation file
variation_lines <- readLines(gzfile(variation_file))

# Function to grep lines that contain any of the gene names
grep_gene_lines <- function(gene_names, lines) {
  matched_lines <- unlist(lapply(gene_names, function(gene) {
    grep(paste0("\\b", gene, "\\b"), lines, value = TRUE)
  }))
  return(unique(matched_lines))  # Ensure unique lines are returned
}

# Get the matched lines
matched_lines <- grep_gene_lines(gene_names, variation_lines)

# Write the matched lines to the output file
writeLines(matched_lines, output_file)

# Print message indicating completion of the first step
cat("Matched lines have been saved to", output_file, "\n")

# Read the matched lines from the output file
matched_lines <- readLines(output_file)

# Filter lines for specific phenotypes
phenotypes <- c("pathogenic", "benign", "uncertain")
grep_phenotype_lines <- function(phenotypes, lines) {
  matched_lines <- unlist(lapply(phenotypes, function(phenotype) {
    grep(phenotype, lines, value = TRUE, ignore.case = TRUE)
  }))
  return(unique(matched_lines))  # Ensure unique lines are returned
}

# Get the lines matching the phenotypes
phenotype_lines <- grep_phenotype_lines(phenotypes, matched_lines)

# Write the phenotype-matched lines to the phenotype output file
writeLines(phenotype_lines, phenotype_output_file)

# Print message indicating completion of the second step
cat("Phenotype-matched lines have been saved to", phenotype_output_file, "\n")
