# Define the file paths as variables
gene_csv_file <- "The_epigenetic_machinery.csv"
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
variation_file <- "homo_sapiens_variation_missense_ClinVar.txt"
output_file <- "calculate_missense_variant_enrichment_within_isoforms_by_pathogenicity_type.txt"

# Load gene names from a CSV file
gene_data <- read.csv(gene_csv_file, header = TRUE)

# Access the Gene_name column values and save to variable
gene_names <- gene_data$Gene_Name

# Append '_HUMAN' to each gene name
gene_names_human <- paste0(gene_names, "_HUMAN")

# Write the modified gene names to a temporary file
writeLines(gene_names_human, "temp_gene_names.txt")

# Read the modified gene names from the temporary file
temp_gene_names <- readLines("temp_gene_names.txt")

# Open the output file for writing
output_conn <- file(output_file, open = "wt")

# Function to safely write lines to the output file
safe_write <- function(line) {
  writeLines(line, output_conn)
}

# Function to count occurrences of a specific variant type, excluding others
count_variants <- function(identifier, variant_type, exclusion_patterns) {
  command <- paste0(
    "grep -i -w '", identifier, "' ", variation_file,
    " | grep -i 'missense variant'",
    " | grep -i '", variant_type, "'",
    " | grep -vi 'somatic'",
    " | grep -vi 'pathogenicity'",
    exclusion_patterns,
    " | wc -l"
  )
  count <- as.numeric(system(command, intern = TRUE))
  return(count)
}

# Define exclusion patterns for each variant type
exclusion_patterns <- list(
  "Likely benign" = " | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "uncertain" = " | grep -vi 'likely benign' | grep -vi 'benign' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Benign" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'likely pathogenic'",
  "Likely pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'pathogenic'"
)

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c ", fasta_file, " | grep '", gene_name_human, "' > temp.txt"))
  
  # Process identifiers and filter out those followed by a hyphen
  system("awk -F '[|]' '{print $2}' temp.txt | awk '!/-/' | sort -u > temp_identifiers.txt")
  temp_identifiers <- readLines("temp_identifiers.txt")
  
  # Loop over each identifier to count mutations for each variant type
  for (identifier in temp_identifiers) {
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      if (count_mutation > 0) {
        safe_write(paste0(identifier, ",", variant_type, ",", count_mutation))
      }
    }
  }
  
  # Extract a different set of identifiers from the same file
  system("awk -F '[|]' '{print $2}' temp.txt | sort -u > temp_identifiers_2.txt")
  temp_identifiers_2 <- readLines("temp_identifiers_2.txt")
  
  # Loop over each identifier to count mutations for each variant type
  for (identifier in temp_identifiers_2) {
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      if (count_mutation > 0) {
        safe_write(paste0(identifier, ",", variant_type, ",", count_mutation))
      }
    }
  }
  
  # Clean up temporary files
  system("rm temp.txt temp_identifiers.txt temp_identifiers_2.txt")
}

# Close the output connection
close(output_conn)

# Remove the temporary gene names file
system("rm temp_gene_names.txt")

# Print completion message
cat("Script execution completed. Results saved to", output_file, "\n")
