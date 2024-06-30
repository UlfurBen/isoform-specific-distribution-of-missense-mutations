# Read the updated amino acid change with property change file
aa_change_data <- readLines("Amino_acid_change_with_property_change_updated.txt")

# Initialize a list to hold the summed counts of property changes
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
  parts <- strsplit(line, " ")[[1]]
  property_change <- parts[2]
  count <- as.numeric(parts[3])
  
  property_change_totals[[property_change]] <- property_change_totals[[property_change]] + count
}

# Prepare the lines for the new file
output_lines <- c()
for (change in names(property_change_totals)) {
  total_count <- property_change_totals[[change]]
  output_lines <- c(output_lines, paste(change, total_count))
}

# Write the summed counts to a new file
writeLines(output_lines, "trial_count_polarity_or_charged_change_property_change_only.txt")

