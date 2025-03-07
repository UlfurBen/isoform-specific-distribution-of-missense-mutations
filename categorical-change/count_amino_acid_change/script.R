# Define file names
input_files <- c(
  "homo_sapiens_variation_missense_ClinVar_benign_likely_benign_only_EM_genes_aa_change.txt",
  "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_EM_genes_aa_change.txt",
  "homo_sapiens_variation_missense_ClinVar_benign_likely_benign_only_all_genes_aa_change.txt",
  "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_aa_change.txt"
)

output_files <- c(
  "homo_sapiens_variation_missense_ClinVar_benign_likely_benign_only_EM_genes_aa_change_count.txt",
  "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_EM_genes_aa_change_count.txt",
  "homo_sapiens_variation_missense_ClinVar_benign_likely_benign_only_all_genes_aa_change_count.txt",
  "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_aa_change_count.txt"
)

# Iterate over each input file and corresponding output file
for (i in 1:length(input_files)) {
  # Read the input file
  data <- read.table(input_files[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get the unique values in the second column and their counts
  aa_change_counts <- table(data[, 2])
  
  # Create a data frame with amino acid change (second column) and their counts
  output_data <- as.data.frame(aa_change_counts)
  
  # Ensure there are 400 rows by filling with missing values if necessary
  all_aa_changes <- unique(data[, 2])
  all_aa_changes <- sort(all_aa_changes)
  
  # Match the original list of amino acid changes with the counts, adding zeros for missing ones
  output_data <- merge(data.frame(Var1 = all_aa_changes), output_data, by = "Var1", all.x = TRUE)
  output_data[is.na(output_data)] <- 0
  
  # Sort by amino acid change (second column)
  output_data <- output_data[order(output_data$Var1), ]
  
  # Write the result to the corresponding output file
  write.table(output_data, output_files[i], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
