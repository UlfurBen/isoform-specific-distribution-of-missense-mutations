# Define the file paths as variables
gene_csv_file <- "The-Epigenetic-Machinery.csv"
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
writeLines("Identifier,Variation count", output_conn)

# Function to safely write lines to the output file
safe_write <- function(line) {
  writeLines(line, output_conn)
}

# Function to count occurrences of missense variants excluding somatic ones
count_variants <- function(identifier, exclude_hyphen = FALSE) {
  exclude_hyphen_pattern <- if (exclude_hyphen) " | grep -v '-'" else ""
  command <- paste0(
    "grep -i -w '", identifier, "' ", variation_file,
    " | grep -i 'missense variant'",
    " | grep -vi 'somatic'",
    exclude_hyphen_pattern,
    " | wc -l"
  )
  print(paste("Running command:", command))  # Debug print statement
  count <- as.numeric(system(command, intern = TRUE))
  print(paste("Count for", identifier, ":", count))  # Debug print statement
  return(count)
}

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  print(paste("Processing gene:", gene_name_human))  # Debug print statement
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c ", fasta_file, " | grep '", gene_name_human, "' > temp.txt"))

  # Process identifiers with hyphens
  system("awk -F '[|]' '{print $2}' temp.txt | grep '-' | sort -u > temp_identifiers_with_hyphen.txt")
  temp_identifiers_with_hyphen <- readLines("temp_identifiers_with_hyphen.txt")
  print(paste("Non-canonical isoforms for", gene_name_human, ":", paste(temp_identifiers_with_hyphen, collapse = ", ")))  # Debug print statement

  # Extract canonical identifiers from non-canonical ones
  canonical_identifiers <- unique(sub("-.*", "", temp_identifiers_with_hyphen))
  print(paste("Canonical isoforms for", gene_name_human, ":", paste(canonical_identifiers, collapse = ", ")))  # Debug print statement

  for (identifier in canonical_identifiers) {
    count_mutation <- count_variants(identifier, exclude_hyphen = TRUE)
    safe_write(paste0(identifier, ",", count_mutation))
  }

  for (identifier in temp_identifiers_with_hyphen) {
    count_mutation <- count_variants(identifier, exclude_hyphen = FALSE)
    safe_write(paste0(identifier, ",", count_mutation))
  }

  # Clean up temporary files
  system("rm temp.txt temp_identifiers_with_hyphen.txt")
}

# Close the output connection
close(output_conn)

# Remove the temporary gene names file
system("rm temp_gene_names.txt")

# Print completion message
cat("Script execution completed. Results saved to", output_file, "\n")
