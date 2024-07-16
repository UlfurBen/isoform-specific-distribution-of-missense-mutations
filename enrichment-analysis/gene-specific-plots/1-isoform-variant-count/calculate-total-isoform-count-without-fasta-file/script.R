# Define the file paths as variables
gene_csv_file <- "The-Epigenetic-Machinery.csv"
variation_file <- "homo_sapiens_variation_missense_ClinVar.txt"
output_file <- "enrichment.txt"

# Load gene names from a CSV file
gene_data <- read.csv(gene_csv_file, header = TRUE)

# Access the Gene_name column values and save to variable
gene_names <- gene_data$Gene_Name

# Open the output file for writing and write the header line
output_conn <- file(output_file, open = "wt")
writeLines("Gene_Name,Identifier,Variation count,Benign,Likely benign,VUS,Likely pathogenic,Pathogenic", output_conn)

# Function to safely write lines to the output file
safe_write <- function(line) {
  writeLines(line, output_conn)
}

# Function to count occurrences of variants excluding somatic ones
count_variants <- function(identifier, exclude_hyphen = FALSE) {
  exclude_hyphen_pattern <- if (exclude_hyphen) " | grep -v '-'" else ""
  command <- paste0(
    "grep -i -w '", identifier, "' ", variation_file,
    " | grep -i 'missense variant'",
    " | grep -vi 'somatic'",
    exclude_hyphen_pattern
  )
  result <- system(command, intern = TRUE)

  # Count specific variant categories
  benign <- length(grep('Benign', result, ignore.case = TRUE))
  likely_benign <- length(grep('Likely benign', result, ignore.case = TRUE))
  vus <- length(grep('uncertain', result, ignore.case = TRUE))
  likely_pathogenic <- length(grep('Likely pathogenic', result, ignore.case = TRUE))
  pathogenic <- length(grep('Pathogenic', result, ignore.case = TRUE))
  
  # Exclude 'pathogenicity'
  pathogenic <- pathogenic - length(grep('pathogenicity', result, ignore.case = TRUE))
  
  total_count <- length(result)
  
  return(list(total_count = total_count, benign = benign, likely_benign = likely_benign, vus = vus, likely_pathogenic = likely_pathogenic, pathogenic = pathogenic))
}

# Process each gene name from the CSV file
for (gene_name in gene_names) {
  print(paste("Processing gene:", gene_name))  # Debug print statement
  
  # Extract identifiers related to the current gene name from the variation file
  command <- paste0("grep -i -w '", gene_name, "' ", variation_file, " | awk -F '\t' '{print $2}' | sort -u > temp_identifiers.txt")
  system(command)
  
  # Read the identifiers from the temporary file
  temp_identifiers <- readLines("temp_identifiers.txt")
  print(paste("Identifiers for", gene_name, ":", paste(temp_identifiers, collapse = ", ")))  # Debug print statement
  
  # Extract canonical identifiers from non-canonical ones
  canonical_identifiers <- unique(sub("-.*", "", temp_identifiers))
  print(paste("Canonical isoforms for", gene_name, ":", paste(canonical_identifiers, collapse = ", ")))  # Debug print statement
  
  for (identifier in canonical_identifiers) {
    variant_counts <- count_variants(identifier, exclude_hyphen = TRUE)
    safe_write(paste0(gene_name, ",", identifier, ",", variant_counts$total_count, ",", variant_counts$benign, ",", variant_counts$likely_benign, ",", variant_counts$vus, ",", variant_counts$likely_pathogenic, ",", variant_counts$pathogenic))
  }
  
  for (identifier in temp_identifiers) {
    variant_counts <- count_variants(identifier, exclude_hyphen = FALSE)
    safe_write(paste0(gene_name, ",", identifier, ",", variant_counts$total_count, ",", variant_counts$benign, ",", variant_counts$likely_benign, ",", variant_counts$vus, ",", variant_counts$likely_pathogenic, ",", variant_counts$pathogenic))
  }
  
  # Clean up temporary files
  system("rm temp_identifiers.txt")
}

# Close the output connection
close(output_conn)

# Print completion message
cat("Script execution completed. Results saved to", output_file, "\n")
