# Load necessary library
library(dplyr)

# Define file paths
input_file <- "enrichment.txt"
output_file <- "enrichment-aggr.txt"

# Read the original enrichment file
data <- read.csv(input_file, header = TRUE)

# Aggregate data by Gene_Name and Identifier
aggregated_data <- data %>%
  group_by(Gene_Name, Identifier) %>%
  summarize(
    Variation_count = sum(Variation.count),
    Benign = sum(Benign),
    Likely_benign = sum(Likely.benign),
    VUS = sum(VUS),
    Likely_pathogenic = sum(Likely.pathogenic),
    Pathogenic = sum(Pathogenic)
  )

# Write the aggregated data to the output file
write.csv(aggregated_data, output_file, row.names = FALSE)

# Print completion message
cat("Aggregation completed. Results saved to", output_file, "\n")
