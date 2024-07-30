# Define the file paths as variables
gene_csv_file <- "The-Epigenetic-Machinery.csv"
fasta_file <- "uniprot_sprot_varsplic.fasta"
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

# Function to aggregate variant counts
aggregate_variants <- function(counts_list) {
  total_count <- sum(sapply(counts_list, function(x) x$total_count))
  benign <- sum(sapply(counts_list, function(x) x$benign))
  likely_benign <- sum(sapply(counts_list, function(x) x$likely_benign))
  vus <- sum(sapply(counts_list, function(x) x$vus))
  likely_pathogenic <- sum(sapply(counts_list, function(x) x$likely_pathogenic))
  pathogenic <- sum(sapply(counts_list, function(x) x$pathogenic))
  
  return(list(total_count = total_count, benign = benign, likely_benign = likely_benign, vus = vus, likely_pathogenic = likely_pathogenic, pathogenic = pathogenic))
}

# Process each modified gene name
for (i in seq_along(temp_gene_names)) {
  gene_name_human <- temp_gene_names[i]
  gene_name <- gene_names[i]
  
  print(paste("Processing gene:", gene_name_human))  # Debug print statement
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("grep '", gene_name_human, "' ", fasta_file, " > temp.txt"))

  # Process identifiers with hyphens
  system("awk -F '[|]' '{print $2}' temp.txt | grep '-' | sort -u > temp_identifiers_with_hyphen.txt")
  temp_identifiers_with_hyphen <- readLines("temp_identifiers_with_hyphen.txt")
  print(paste("Non-canonical isoforms for", gene_name_human, ":", paste(temp_identifiers_with_hyphen, collapse = ", ")))  # Debug print statement

  # Extract canonical identifiers from non-canonical ones
  canonical_identifiers <- unique(sub("-.*", "", temp_identifiers_with_hyphen))
  print(paste("Canonical isoforms for", gene_name_human, ":", paste(canonical_identifiers, collapse = ", ")))  # Debug print statement

  # Collect counts for canonical identifiers
  canonical_counts <- list()
  for (identifier in canonical_identifiers) {
    variant_counts <- count_variants(identifier, exclude_hyphen = TRUE)
    canonical_counts[[identifier]] <- variant_counts
  }
  
  # Collect counts for non-canonical identifiers
  non_canonical_counts <- list()
  for (identifier in temp_identifiers_with_hyphen) {
    variant_counts <- count_variants(identifier, exclude_hyphen = FALSE)
    non_canonical_counts[[identifier]] <- variant_counts
  }

  # Aggregate counts for canonical identifiers
  for (identifier in canonical_identifiers) {
    if (identifier %in% names(non_canonical_counts)) {
      counts_list <- list(canonical_counts[[identifier]], non_canonical_counts[[identifier]])
      aggregated_counts <- aggregate_variants(counts_list)
      safe_write(paste0(gene_name, ",", identifier, ",", aggregated_counts$total_count, ",", aggregated_counts$benign, ",", aggregated_counts$likely_benign, ",", aggregated_counts$vus, ",", aggregated_counts$likely_pathogenic, ",", aggregated_counts$pathogenic))
    } else {
      counts <- canonical_counts[[identifier]]
      safe_write(paste0(gene_name, ",", identifier, ",", counts$total_count, ",", counts$benign, ",", counts$likely_benign, ",", counts$vus, ",", counts$likely_pathogenic, ",", counts$pathogenic))
    }
  }

  # Write counts for non-canonical identifiers not in canonical list
  for (identifier in temp_identifiers_with_hyphen) {
    if (!(identifier %in% canonical_identifiers)) {
      counts <- non_canonical_counts[[identifier]]
      safe_write(paste0(gene_name, ",", identifier, ",", counts$total_count, ",", counts$benign, ",", counts$likely_benign, ",", counts$vus, ",", counts$likely_pathogenic, ",", counts$pathogenic))
    }
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
