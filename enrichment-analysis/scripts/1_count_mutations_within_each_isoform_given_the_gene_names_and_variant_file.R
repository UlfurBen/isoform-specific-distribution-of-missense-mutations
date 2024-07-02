# Load gene names from a CSV file
gene_data <- read.csv("The_epigenetic_machinery.csv", header = TRUE)

# Output file
output_file <- "calculate_missense_variant_enrichment_within_isoforms.txt"

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

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract gene uniprot identifiers from the FASTA file
  system(paste0("gunzip -c uniprot_sprot_varsplic.fasta.gz | grep '", gene_name_human, "' > temp.txt"))
  
  # Process identifiers
  system("awk -F '[|-]' '{print $2}' temp.txt | sort -u > temp_identifiers.txt")
  temp_identifiers <- readLines("temp_identifiers.txt")
  
  # Loop over each identifier to count mutations, only taking the first isoform which specifically is not followed by a hyphen
  for (identifier in temp_identifiers) {
    count_mutation <- as.numeric(system(paste0("awk '$0 ~ /", identifier, "/ && $0 !~ /", identifier, "-/ && $0 ~ /missense variant/ {c++} END {print c+0}' filter_database_by_EM_genes_missense_variants_ClinVar_with_RCV.txt"), intern = TRUE))
    if (count_mutation > 0) {
      write(paste0(identifier, ",", count_mutation, "\n"), file = output_conn, append = TRUE)
    }
  }
  
  # Extract a different set of identifiers from the same file
  system("awk -F '|' '{print $2}' temp.txt | sort -u > temp_identifiers_2.txt")
  temp_identifiers_2 <- readLines("temp_identifiers_2.txt")
  
  # Loop over each identifier to count mutations, specifically taking only the second to last isoform, i.e. isoforms with a hyphen
  for (identifier in temp_identifiers_2) {
    count_mutation <- as.numeric(system(paste0("grep -w '", identifier, "' filter_database_by_EM_genes_missense_variants_ClinVar_with_RCV.txt | grep 'missense variant' | wc -l"), intern = TRUE))
    if (count_mutation > 0) {
      write(paste0(identifier, ",", count_mutation, "\n"), file = output_conn, append = TRUE)
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
