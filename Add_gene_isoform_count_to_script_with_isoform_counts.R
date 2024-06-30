# Define file paths
input_file <- "merged_calculate_missense_variant_enrichment_within_isoforms_with_lengths_all_300.txt"
# On elja -> ~/uniprot_sprot_varsplic.fasta.gz
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
output_file <- "appended_with_isoform_counts.txt"


# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Function to count the number of isoforms for a given gene from the FASTA file
count_isoforms <- function(gene_id, fasta_file) {
  # Read the FASTA file
  fasta_lines <- readLines(gzfile(fasta_file))

  # Initialize the isoform count
  isoform_count <- 1

  # Loop through each line in the FASTA file
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      # Check if the line contains the gene_id
      if (grepl(paste0("\\|", gene_id, "-"), line)) {
        isoform_count <- isoform_count + 1
      }
    }
  }

  return(isoform_count)
}

# Extract the unique gene IDs from the data
gene_ids <- unique(sub("-.*", "", data$Isoform))

# Create a list to store the isoform counts
isoform_counts <- sapply(gene_ids, count_isoforms, fasta_file = fasta_file)

# Append the isoform count to each line in the data
data$Isoform_Count <- sapply(sub("-.*", "", data$Isoform), function(x) isoform_counts[x])

# Write the updated data to the output file
write.table(data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print a message indicating completion
cat("Data with appended isoform counts has been saved to", output_file, "\n")
