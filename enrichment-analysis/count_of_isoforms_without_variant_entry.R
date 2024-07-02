# Read the list of distinct isoform IDs
isoform_ids <- readLines("list_of_distinct_isoform_ids.txt")

# Read the homo sapiens variation missense ClinVar data
clinvar_data <- readLines("homo_sapiens_variation_missense_ClinVar_Reference_SNP.txt")

# Convert ClinVar data into a single string for faster searching
clinvar_data_combined <- paste(clinvar_data, collapse = " ")

# Initialize a counter for entries without matches
no_match_counter <- 0

# Create a list to store the count messages
count_messages <- c()

# Loop over each isoform ID and check for matches in the ClinVar data
for (isoform_id in isoform_ids) {
  # Check if the isoform ID is found in the combined ClinVar data string
  if (!grepl(isoform_id, clinvar_data_combined)) {
    no_match_counter <- no_match_counter + 1

    # Every 100 increments, store the count message
    if (no_match_counter %% 100 == 0) {
      count_messages <- c(count_messages, paste0("there are: ", no_match_counter))
    }
  }
}

# Store the final count if it's not a multiple of 100
if (no_match_counter %% 100 != 0) {
  count_messages <- c(count_messages, paste0("there are: ", no_match_counter))
}

# Write all count messages to the file in one operation
write(count_messages, file = "count_of_isoforms_without_variant_entry_quicket.txt")

# Print the number of entries in the first file that don't have matches in the second file
cat("Number of entries in the first file without matches in the second file:", no_match_counter, "\n")

