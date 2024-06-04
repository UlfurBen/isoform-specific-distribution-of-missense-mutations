# Load required libraries
library(ggplot2)
library(Biostrings)
library(dplyr)
library(stringr)

# Load mutation data from "kmt2d-pathogenic-mutations.txt"
mutations <- read.table("kmt2d-pathogenic-mutations.txt", header = TRUE, sep = "\t")

# Convert the data to a data frame
mutations_df <- data.frame(mutations)

# Extract the third column as a list
extract_protein_change <- function(df) {
  # Extract the third column
  protein_changes <- df[, 3]
  
  # Convert to list
  result <- as.list(protein_changes)
  
  return(result)
}

# Extracted list of protein changes
protein_changes_list <- extract_protein_change(mutations_df)

# Print the result
print(protein_changes_list)

# Save the protein_changes_list to a text file
writeLines(unlist(protein_changes_list), "protein_change.txt")


print(mutations_df)

# Convert the "Name" column to a factor
mutations_df$Name <- factor(mutations_df$Name)

# Load isoform data from "AIRE_isoform_sequence-1.fasta"
isoform_files <- c("AIRE_isoform_sequence-1.fasta")
print(isoform_files)
readDNAStringSet(file="AIRE_isoform_sequence-1.fasta", format="fasta")

isoform_seqs <- DNAStringSet()

for (file in isoform_files) {
  isoform_seqs <- c(isoform_seqs, readDNAStringSet(file))
}


