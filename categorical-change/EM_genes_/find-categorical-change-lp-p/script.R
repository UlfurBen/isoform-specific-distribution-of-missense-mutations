# Define file paths
# Get input_file by running grep -i 'pathogenic' homo_sapiens_variation_missense_ClinVar.txt | grep -i 'likely pathogenic' | grep -vi 'uncertain' | grep -vi 'benign' | grep -vi 'likely benign' > homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only.txt
input_file <- "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only.txt"
property_lookup_file <- "Amino_acid_change_with_property_change.txt"
identifier_file <- "The-Epigenetic-Machinery.csv"

# Output
output_file <- "categorical-change-EM-genes.txt"

# Read the input files
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
properties <- read.table(property_lookup_file, header = FALSE, sep = " ", stringsAsFactors = FALSE, quote = "")
identifiers <- read.csv(identifier_file, header = TRUE, stringsAsFactors = FALSE)

# Rename columns for better understanding
colnames(data) <- c("Column1", "Isoform", "AminoAcidChange", "Column4", "Column5", "Column6", "Column7", "Column8", "Column9", "Column10", "Column11", "Column12", "Identifier", "Column14")
colnames(properties) <- c("AminoAcidChange", "PropertyChange")

# Extract relevant identifiers from the identifier file
relevant_identifiers <- identifiers$Identifier

# Filter the data for relevant identifiers in the 13th column
filtered_data <- data[data$Identifier %in% relevant_identifiers, ]

# Function to remove numeric characters between amino acid names
filter_amino_acids <- function(variant) {
  gsub("p\\.([A-Za-z]+)\\d+([A-Za-z]+)", "\\1\\2", variant)
}

# Apply the function to the relevant column (assuming the third column contains the amino acid changes)
filtered_data$AminoAcidChange <- sapply(filtered_data$AminoAcidChange, filter_amino_acids)

# Merge the data with properties based on the amino acid change
merged_data <- merge(filtered_data, properties, by = "AminoAcidChange")

# Select only the Isoform, AminoAcidChange, and PropertyChange columns for output
output_data <- merged_data[, c("Isoform", "AminoAcidChange", "PropertyChange")]

# Write the modified data to the output file
write.table(output_data, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Print a message indicating completion
cat("Filtered data with property changes has been saved to", output_file, "\n")