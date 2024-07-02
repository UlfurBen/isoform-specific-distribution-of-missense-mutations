# Input 
input <- variation_amino_acid_change_with_property_change.txt
# Read the filtered variant file
filtered_variants <- readLines(input)

# Output
output <- Amino_acid_change_with_property_change_updated.txt

# Read in the list of all possible amino acid changes with their correspongding property changes
aa_change_properties <- readLines("Amino_acid_change_with_property_change.txt")

# Create a map for property changes with counts initialized to 0
property_change_counts <- list()
for (line in aa_change_properties) {
  parts <- strsplit(line, " ")[[1]]
  aa_pair <- parts[1]
  property_change <- parts[2]
  property_change_counts[[aa_pair]] <- 0
}

# Count the occurrences of each amino acid change in the filtered variants
for (line in filtered_variants) {
  fields <- strsplit(line, " ")[[1]]
  aa_change <- fields[2]
  
  if (aa_change %in% names(property_change_counts)) {
    property_change_counts[[aa_change]] <- property_change_counts[[aa_change]] + 1
  }
}

# Prepare the lines to append to the Amino_acid_change_with_property_change.txt file
output_lines <- c()
for (line in aa_change_properties) {
  parts <- strsplit(line, " ")[[1]]
  aa_pair <- parts[1]
  property_change <- parts[2]
  count <- property_change_counts[[aa_pair]]
  output_lines <- c(output_lines, paste(aa_pair, property_change, count))
}

# Write the updated counts to a new file
writeLines(output_lines, output)

