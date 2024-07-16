# Input
input <- 'homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_aa_change.txt'

# Output
output <- 'homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_categorical-count.txt'

# Read the updated amino acid change with property change file
aa_change_data <- readLines(input)

# Initialize a list to hold the counts of property changes
property_change_totals <- list(
  "Nonpolar->Nonpolar" = 0,
  "Nonpolar->Polar" = 0,
  "Nonpolar->Negatively_charged" = 0,
  "Nonpolar->Positively_charged" = 0,
  "Polar->Nonpolar" = 0,
  "Polar->Polar" = 0,
  "Polar->Negatively_charged" = 0,
  "Polar->Positively_charged" = 0,
  "Negatively_charged->Nonpolar" = 0,
  "Negatively_charged->Polar" = 0,
  "Negatively_charged->Negatively_charged" = 0,
  "Negatively_charged->Positively_charged" = 0,
  "Positively_charged->Nonpolar" = 0,
  "Positively_charged->Polar" = 0,
  "Positively_charged->Negatively_charged" = 0,
  "Positively_charged->Positively_charged" = 0
)

# Sum the counts for each property change type
for (line in aa_change_data) {
  parts <- strsplit(line, "\\s+")[[1]]
  if (length(parts) >= 3) {
    property_change <- parts[3]
    
    if (!is.null(property_change_totals[[property_change]])) {
      property_change_totals[[property_change]] <- property_change_totals[[property_change]] + 1
    } else {
      warning(paste("Invalid property change in line:", line))
    }
  } else {
    warning(paste("Line does not have enough parts:", line))
  }
}

# Prepare the lines for the new file
output_lines <- c()
for (change in names(property_change_totals)) {
  total_count <- property_change_totals[[change]]
  output_lines <- c(output_lines, paste(change, total_count))
}

# Write the summed counts to a new file
writeLines(output_lines, output)
