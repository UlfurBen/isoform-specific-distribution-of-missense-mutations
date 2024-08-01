# Load necessary library
library(dplyr)

# Read the BED file
bed_data <- read.table("filtered_no_scientific_notation.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Filter rows where any column contains 'crebbp' (case insensitive)
filtered_data <- bed_data %>% 
  filter_all(any_vars(grepl("crebbp", ., ignore.case = TRUE)))

# Extract unique values from the 8th column
unique_values <- unique(filtered_data$V8)

# Write unique values to a new BED file
write.table(unique_values, file = "crebbp_unique_ENST.bed", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create the main folder 'crebbp'
if (!dir.exists("crebbp")) {
  dir.create("crebbp")
}

# Create subfolders for each unique ENST value
for (value in unique_values) {
  dir.create(file.path("crebbp", value))
}

# Output message to confirm completion
cat("Folders created successfully.\n")
