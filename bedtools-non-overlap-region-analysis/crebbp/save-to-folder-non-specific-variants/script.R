# Load necessary libraries
library(dplyr)

# Read the BED file
bed_data <- read.table("homo_sapiens_variation_missense_ClinVar_filtered_relevancy_pathogenic.bed
", header = FALSE, stringsAsFactors = FALSE, fill=TRUE)

# Filter the data for 'crebbp' (case insensitive) in the fifth column
filtered_data <- bed_data %>%
  filter(grepl("crebbp", V5, ignore.case = TRUE))

# Create the main folder
main_folder <- "crebbp"
if (!dir.exists(main_folder)) {
  dir.create(main_folder)
}

# Find all unique values in the entire row that start with "ENST"
unique_values <- unique(unlist(filtered_data %>% select(everything()) %>% as.vector()))
unique_values <- unique_values[grepl("^ENST", unique_values)]

# Save the filtered rows in the main folder
for (value in unique_values) {
  subset_data <- filtered_data %>% filter_all(any_vars(grepl(value, .)))
  write.table(subset_data, file = file.path(main_folder, paste0(value, ".bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Print a completion message
cat("Processing complete. Filtered data saved in", main_folder, "\n")
