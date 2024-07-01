# Define the file paths
fasta_file <- "uniprot_sprot_varsplic.fasta.gz"
variant_file <- "calculate_missense_variant_enrichment_within_isoforms.txt"
output_file <- "calculate_missense_variant_enrichment_within_isoforms_with_lengths.txt"

# Read the entire fasta file into memory
fasta_lines <- readLines(gzfile(fasta_file))

# Read the variant file
variant_data <- read.csv(variant_file, header = FALSE)
colnames(variant_data) <- c("Isoform", "Count")

# Initialize an empty list to store isoform sequences
isoform_sequences <- list()

# Initialize variables to store the current isoform id and sequence
current_isoform <- NULL
current_sequence <- ""

# Loop through the fasta lines to extract isoform sequences
for (line in fasta_lines) {
  if (startsWith(line, ">")) {
    # If there's an ongoing isoform, save its sequence
    if (!is.null(current_isoform)) {
      isoform_sequences[[current_isoform]] <- current_sequence
    }
    # Extract the isoform id from the header line
    current_isoform <- sub(".*\\|(.*?)\\|.*", "\\1", line)
    current_sequence <- ""
  } else {
    # Append the sequence line
    current_sequence <- paste0(current_sequence, line)
  }
}
# Save the last isoform sequence
if (!is.null(current_isoform)) {
  isoform_sequences[[current_isoform]] <- current_sequence
}

# Function to calculate sequence length
calculate_length <- function(sequence) {
  return(nchar(sequence))
}

# Calculate lengths for the isoforms in the variant file
variant_data$Length <- sapply(variant_data$Isoform, function(isoform) {
  if (isoform %in% names(isoform_sequences)) {
    return(calculate_length(isoform_sequences[[isoform]]))
  } else {
    return(NA)
  }
})

# Write the results to the output file
write.table(variant_data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print message indicating completion
cat("Isoform lengths have been appended to", output_file, "\n")

