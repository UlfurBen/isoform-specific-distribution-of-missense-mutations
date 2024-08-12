# Load necessary library
library(dplyr)

# Define file paths
format_results_path <- "~/final-results/format_results.bed"
benign_format_results_path <- "~/final-results/benign_format_results.bed"

# Read the files without headers
format_results <- read.table(format_results_path, header = FALSE, sep = "\t")
benign_format_results <- read.table(benign_format_results_path, header = FALSE, sep = "\t")

# Perform the join based on the first 3 columns
merged_results <- left_join(format_results, benign_format_results[, c(1:3, 6)], 
                            by = c("V1" = "V1", "V2" = "V2", "V3" = "V3"))

# Write the result to a new file called appended-results.bed
output_path <- "~/final-results/appended-results.bed"
write.table(merged_results, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Optional: Print a message indicating completion
cat("Merging complete. Results saved to", output_path, "\n")
