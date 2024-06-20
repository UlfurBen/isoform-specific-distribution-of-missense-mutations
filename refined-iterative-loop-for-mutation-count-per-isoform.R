# Load gene names from a CSV file
gene_data <- read.csv("The-Epigenetic-Machinery.csv", header = TRUE)

# Assuming the gene names are in a column named 'Gene_Name'
gene_names <- gene_data$Gene_Name

# Append '_HUMAN' to each gene name
gene_names_human <- paste0(gene_names, "_HUMAN")

# Write the modified gene names to a temporary file
writeLines(gene_names_human, "temp_gene_names.txt")

# Read the modified gene names from the temporary file
temp_gene_names <- readLines("temp_gene_names.txt")
output_file <- "output.txt"

# Process each modified gene name
for (gene_name_human in temp_gene_names) {
  # Extract identifiers from the FASTA file
  system(paste0("gunzip -c uniprot_sprot_varsplic.fasta.gz | grep '", gene_name_human, "' > temp.txt"))

# Process identifiers and count occurrences
system("awk -F '[|-]' '{print $2}' temp.txt | sort -u > temp_identifiers.txt")
temp_identifiers <- readLines("temp_identifiers.txt")

# Loop over each identifier to count mutations and gnomAD entries
for (identifier in temp_identifiers) {
  count_mutation <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation.txt.gz | awk '$0 ~ /", identifier, "/ && $0 !~ /", identifier, "-/ && $0 ~ /missense variant/ {c++} END {print c+0}'"), intern = TRUE))
  count_gnomAD <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation.txt.gz | awk '$0 ~ /", identifier, "/ && $0 !~ /", identifier, "-/ && $0 ~ /gnomAD/ && $0 ~ /missense variant/ {c++} END {print c+0}'"), intern = TRUE))
  write(paste0(identifier, ",", count_mutation, ",", count_gnomAD), file = output_file, append = TRUE)
}
    # Extract a different set of identifiers from the same file
  system("awk -F '|' '{print $2}' temp.txt | sort -u > temp_identifiers_2.txt")
  temp_identifiers_2 <- readLines("temp_identifiers_2.txt")

  # Loop over each identifier to count mutations and gnomAD entries
  for (identifier in temp_identifiers_2) {
    count_mutation <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation.txt.gz | grep -w '", identifier, "' | grep 'missense variant' | wc -l"), intern = TRUE))
    count_gnomAD <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation.txt.gz | grep -w '", identifier, "' | grep 'gnomAD' | grep 'missense variant' | wc -l"), intern = TRUE))
    write(paste0(identifier, ",", count_mutation, ",", count_gnomAD), file = output_file, append = TRUE)
  }

  # Clean up temporary files
  system("rm temp.txt temp_identifiers.txt temp_identifiers_2.txt")
}

# Remove the temporary gene names file
system("rm temp_gene_names.txt")
