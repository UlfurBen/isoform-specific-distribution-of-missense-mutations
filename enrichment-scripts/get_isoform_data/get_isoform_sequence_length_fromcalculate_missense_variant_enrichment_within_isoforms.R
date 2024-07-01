# Define file paths
variants_file <- "calculate_missense_variant_enrichment_within_isoforms.txt"
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
final_output_file <- "merged_isoform_lengths_with_mutation_count.txt"

# Read the variants file
variants_data <- read.csv(variants_file, header = FALSE, stringsAsFactors = FALSE)
colnames(variants_data) <- c("Isoform", "Count")

# Function to read gzipped FASTA file manually
read_fasta_gz <- function(gz_file) {
  lines <- readLines(gzfile(gz_file))
  fasta_sequences <- list()
  seq_name <- NULL
  seq <- NULL
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(seq_name)) {
        fasta_sequences[[seq_name]] <- paste(seq, collapse = "")
      }
      seq_name <- sub(" .*", "", sub("^>", "", line))
      seq <- c()
    } else {
      seq <- c(seq, line)
    }
  }
  if (!is.null(seq_name)) {
    fasta_sequences[[seq_name]] <- paste(seq, collapse = "")
  }
  return(fasta_sequences)
}

# Function to get sequence length from UniProt REST API
get_sequence_length_from_uniprot <- function(isoform_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", isoform_id, ".fasta")
  fasta_data <- readLines(url)
  seq <- paste(fasta_data[!startsWith(fasta_data, ">")], collapse = "")
  return(nchar(seq))
}

# Read the FASTA file
fasta_sequences <- read_fasta_gz(fasta_file)

# Initialize lists to store lengths for isoforms
isoform_lengths <- list()

# Process each isoform in the variants data
for (i in 1:nrow(variants_data)) {
  isoform_id <- variants_data$Isoform[i]
  length <- NA
  
  # Check if the isoform is in the FASTA file
  if (isoform_id %in% names(fasta_sequences)) {
    length <- nchar(fasta_sequences[[isoform_id]])
  } else {
    # For isoforms without a hyphen, fetch from UniProt REST API
    if (!grepl("-", isoform_id)) {
      length <- get_sequence_length_from_uniprot(isoform_id)
    }
  }
  
  # Store the length in the list
  isoform_lengths[[isoform_id]] <- length
}

# Add the lengths to the variants data
variants_data$Length <- sapply(variants_data$Isoform, function(id) isoform_lengths[[id]])

# Filter out rows with NA values in the Length column
filtered_variants_data <- na.omit(variants_data)

# Write the merged results to the final output file
write.table(filtered_variants_data, file = final_output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print message indicating completion
cat("Isoform lengths have been processed, merged, and saved to", final_output_file, "\n")
