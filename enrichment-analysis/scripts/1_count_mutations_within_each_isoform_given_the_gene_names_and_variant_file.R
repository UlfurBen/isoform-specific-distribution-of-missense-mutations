# Define the file paths as variables
gene_csv_file <- "The_Epigenetic_Machinery.csv"
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
variation_file <- "homo_sapiens_variation_missense_ClinVar.txt"
output_file <- "enrichment.txt"

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

# Open the output file for writing and write the header line
output_conn <- file(output_file, open = "wt")
writeLines("Identifier,Likely benign,uncertain,Benign,Pathogenic,Likely pathogenic", output_conn)

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
  "Likely benign" = " | grep -vi 'uncertain' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "uncertain" = " | grep -vi 'likely benign' | grep -vi 'benign' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Benign" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'likely pathogenic'",
  "Likely pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign'"
)

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c ", fasta_file, " | grep '", gene_name_human, "' > temp.txt"))
  
  # Process identifiers without hyphens (canonical isoforms)
  system("awk -F '[|]' '{print $2}' temp.txt | grep -v '-' | sort -u > temp_identifiers_no_hyphen.txt")
  temp_identifiers_no_hyphen <- readLines("temp_identifiers_no_hyphen.txt")
  
  for (identifier in temp_identifiers_no_hyphen) {
    variant_counts <- c()
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      variant_counts <- c(variant_counts, count_mutation)
    }
    safe_write(paste0(identifier, ",", paste(variant_counts, collapse = ",")))
  }

  # Process identifiers with hyphens
  system("awk -F '[|]' '{print $2}' temp.txt | grep '-' | sort -u > temp_identifiers_with_hyphen.txt")
  temp_identifiers_with_hyphen <- readLines("temp_identifiers_with_hyphen.txt")
  
  for (identifier in temp_identifiers_with_hyphen) {
    variant_counts <- c()
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      variant_counts <- c(variant_counts, count_mutation)
    }
    safe_write(paste0(identifier, ",", paste(variant_counts, collapse = ",")))
  }
  
  # Clean up temporary files
  system("rm temp.txt temp_identifiers_no_hyphen.txt temp_identifiers_with_hyphen.txt")
}

# Close the output connection
close(output_conn)

# Remove the temporary gene names file
system("rm temp_gene_names.txt")

# Print completion message
cat("Script execution completed. Results saved to", output_file, "\n")


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

# Open the output file for writing and write the header line
output_conn <- file(output_file, open = "wt")
writeLines("Identifier,Likely benign,uncertain,Benign,Pathogenic,Likely pathogenic", output_conn)

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
  "Likely benign" = " | grep -vi 'uncertain' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "uncertain" = " | grep -vi 'likely benign' | grep -vi 'benign' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Benign" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'pathogenic' | grep -vi 'likely pathogenic'",
  "Pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'likely pathogenic'",
  "Likely pathogenic" = " | grep -vi 'likely benign' | grep -vi 'uncertain' | grep -vi 'benign'"
)

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c ", fasta_file, " | grep '", gene_name_human, "' > temp.txt"))
  
  # Process identifiers without hyphens (canonical isoforms)
  system("awk -F '[|]' '{print $2}' temp.txt | grep -v '-' | sort -u > temp_identifiers_no_hyphen.txt")
  temp_identifiers_no_hyphen <- readLines("temp_identifiers_no_hyphen.txt")
  
  for (identifier in temp_identifiers_no_hyphen) {
    variant_counts <- c()
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      variant_counts <- c(variant_counts, count_mutation)
    }
    safe_write(paste0(identifier, ",", paste(variant_counts, collapse = ",")))
  }

  # Process identifiers with hyphens
  system("awk -F '[|]' '{print $2}' temp.txt | grep '-' | sort -u > temp_identifiers_with_hyphen.txt")
  temp_identifiers_with_hyphen <- readLines("temp_identifiers_with_hyphen.txt")
  
  for (identifier in temp_identifiers_with_hyphen) {
    variant_counts <- c()
    for (variant_type in names(exclusion_patterns)) {
      count_mutation <- count_variants(identifier, variant_type, exclusion_patterns[[variant_type]])
      variant_counts <- c(variant_counts, count_mutation)
    }
    safe_write(paste0(identifier, ",", paste(variant_counts, collapse = ",")))
  }
  
  # Clean up temporary files
  system("rm temp.txt temp_identifiers_no_hyphen.txt temp_identifiers_with_hyphen.txt")
}

# Close the output connection
close(output_conn)

# Remove the temporary gene names file
system("rm temp_gene_names.txt")

# Print completion message
cat("Script execution completed. Results saved to", output_file, "\n")
