# Define the input file path
input_file <- "uniprot_sprot_varsplic.fasta.gz"

# Use system command to gunzip and read the contents of the file
fasta_content <- system(paste("gunzip -c", input_file), intern = TRUE)

# Initialize a vector to store isoform IDs
isoform_ids <- c()

# Loop through each line in the fasta_content
for (line in fasta_content) {
  if (startsWith(line, ">")) {
    # Extract the isoform ID using regular expression
    match <- regmatches(line, regexec(">sp\\|([^|]+)\\|", line))
    if (length(match[[1]]) > 1) {
      isoform_id <- match[[1]][2]
      isoform_ids <- c(isoform_ids, isoform_id)
    }
  }
}

# Get the unique isoform IDs
unique_isoform_ids <- unique(isoform_ids)

# Write the unique isoform IDs to a file
output_file <- "list_of_distinct_isoform_ids.txt"
write(unique_isoform_ids, file = output_file)

# Print the number of unique isoform IDs
print(paste("Number of unique isoform IDs:", length(unique_isoform_ids)))

