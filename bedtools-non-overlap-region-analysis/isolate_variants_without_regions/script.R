# Load necessary library
library(dplyr)

# Read the BED file
bed_data <- read.table("homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Filter rows where the 4th column contains 'crebbp' (case insensitive)
filtered_data <- bed_data %>%
  filter(grepl("crebbp", V4, ignore.case = TRUE))

# Extract unique values from the 8th column
unique_values <- unique(filtered_data$V8)

# Create the main folder 'crebbp'
if (!dir.exists("crebbp")) {
  dir.create("crebbp")
}

# Create subfolders for each unique ENST value and save respective rows
for (value in unique_values) {
  subfolder_path <- file.path("crebbp", value)
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path)
  }
  
  # Filter rows that correspond to the current unique value in the 8th column
  subset_data <- filtered_data %>%
    filter(V8 == value)
  
  # Write the filtered data to a file in the respective subfolder
  write.table(subset_data, file = file.path(subfolder_path, paste0(value, ".bed")), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

# Output message to confirm completion
cat("Folders and files created successfully.\n")
