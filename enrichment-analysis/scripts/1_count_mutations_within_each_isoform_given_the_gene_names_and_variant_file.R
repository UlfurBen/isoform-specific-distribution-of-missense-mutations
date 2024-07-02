# Define the file paths as variables
gene_csv_file <- "The_epigenetic_machinery.csv"
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
variation_file <- "homo_sapiens_variation_missense_ClinVar.txt"
output_file <- "calculate_pathogenic_or_likely_pathogenic_variants_per_isoform.txt"

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

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c ", fasta_file, " | grep '", gene_name_human, "' > temp.txt"))
  
  # Process identifiers
  system("awk -F '[|-]' '{print $2}' temp.txt | sort -u > temp_identifiers.txt")
  temp_identifiers <- readLines("temp_identifiers.txt")
  
  # Loop over each identifier to count mutations, only taking the first isoform which specifically is not followed by a hyphen
  for (identifier in temp_identifiers) {
    count_mutation <- as.numeric(system(paste0("awk '$0 ~ /", identifier, "/ && $0 !~ /", identifier, "-/ && $0 ~ /missense variant/ && ($0 ~ /Pathogenic/ || $0 ~ /Likely pathogenic/) && $0 !~ /pathogenicity/ && $0 !~ /uncertain/ && $0 !~ /benign/ {c++} END {print c+0}' ", variation_file), intern = TRUE))
    if (count_mutation > 0) {
      safe_write(paste0(identifier, ",", count_mutation))
    }
  }
  
  # Extract a different set of identifiers from the same file
  system("awk -F '|' '{print $2}' temp.txt | sort -u > temp_identifiers_2.txt")
  temp_identifiers_2 <- readLines("temp_identifiers_2.txt")
  
  # Loop over each identifier to count mutations, specifically taking only the second to last isoform, i.e. isoforms with a hyphen
  for (identifier in temp_identifiers_2) {
    count_mutation <- as.numeric(system(paste0("grep -w '", identifier, "' ", variation_file, " | grep 'missense variant' | grep -E 'Pathogenic|Likely pathogenic' | grep -v 'pathogenicity' | grep -v 'uncertain' | grep -v 'benign' | wc -l"), intern = TRUE))
    if (count_mutation > 0) {
      safe_write(paste0(identifier, ",", count_mutation))
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
