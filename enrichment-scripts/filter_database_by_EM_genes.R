# Define the file paths
gene_csv_file <- "The-Epigenetic-Machinery.csv"
clinvar_file <- "homo_sapiens_variation_missense_ClinVar_Reference_SNP.txt"
output_file <- "homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt"

# Read gene names from the TXT file which can be handled like a CSV file
gene_data <- read.csv(gene_csv_file, header = TRUE)
gene_names <- gene_data$Gene_Name

# Read the lines of the ClinVar SNP file
clinvar_lines <- readLines(clinvar_file)

# Function to grep lines that contain any of the gene names
grep_gene_lines <- function(gene_names, lines) {
  matched_lines <- unlist(lapply(gene_names, function(gene) {
    grep(paste0("\\b", gene, "\\b"), lines, value = TRUE)
  }))
  return(unique(matched_lines))  # Ensure unique lines are returned
}

# Get the matched lines
matched_lines <- grep_gene_lines(gene_names, clinvar_lines)

# Write the matched lines to the output file
writeLines(matched_lines, output_file)

# Print message indicating completion
cat("Matched lines have been saved to", output_file, "\n")

