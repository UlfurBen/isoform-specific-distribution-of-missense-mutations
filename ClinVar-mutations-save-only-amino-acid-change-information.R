# Load required libraries
library(ggplot2)
library(Biostrings)
library(dplyr)
library(stringr)

# Load mutation data from "kmt2d-pathogenic-mutations.txt"
mutations <- read.table("kmt2d-pathogenic-mutations.txt", header = TRUE, sep = "\t")

# Convert the data to a data frame
mutations_df <- data.frame(mutations)

# Function to extract amino acid change information and save to a list
extract_mutations <- function(df) {
  # Regular expression to match the "c.16442G>A" part
  pattern <- "c\\.[0-9]+[A-Z]>[A-Z]"
  
  # Extract matches
  matches <- str_extract(df[, 1], pattern)
  
  # Convert matches to list
  result <- as.list(matches)
  
  return(result)
}

# Extracted list of mutations
mutations_list <- extract_mutations(mutations_df)

# Print the result
print(mutations_list)
writeLines(unlist(mutations_list), "mutations_list.txt")
