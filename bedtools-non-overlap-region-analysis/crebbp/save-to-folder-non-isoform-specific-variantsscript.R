# Load necessary libraries
library(dplyr)

# Read the BED file
bed_data <- read.table("filtered_no_scientific_notation.bed", header = FALSE, stringsAsFactors = FALSE)

# Filter the data for 'crebbp' (case insensitive) in the fifth column
filtered_data <- bed_data %>%
  filter(grepl("crebbp", V5, ignore.case = TRUE))

# Create the main folder
main_folder <- "crebbp"
if (!dir.exists(main_folder)) {
  dir.create(main_folder)
}

# Get unique values from the fourth column
unique_values <- unique(filtered_data$V4)

# Save the filtered rows in the main folder
for (value in unique_values) {
  subset_data <- filtered_data %>% filter(V4 == value)
  write.table(subset_data, file = file.path(main_folder, paste0(value, ".bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Print a completion message
cat("Processing complete. Filtered data saved in", main_folder, "\n")
