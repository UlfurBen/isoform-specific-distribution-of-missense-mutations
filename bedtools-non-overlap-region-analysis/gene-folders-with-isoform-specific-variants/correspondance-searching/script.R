# After running the script run: column -t ~/final-results/appended-results.bed > ~/final-results/appended-formatted-results.bed

# Load necessary library
library(dplyr)

# Define file paths
results_path <- "~/final-results/results.bed"
benign_results_path <- "~/final-results/benign_results.bed"

# Read the files without headers
results <- read.table(results_path, header = FALSE, sep = "\t")
benign_results <- read.table(benign_results_path, header = FALSE, sep = "\t")

# Perform the join based on the first 3 columns
merged_results <- left_join(results, benign_results[, c(1:3, 6)], 
                            by = c("V1" = "V1", "V2" = "V2", "V3" = "V3"))

# Write the result to a new file called appended-results.bed
output_path <- "~/final-results/appended-results.bed"
write.table(merged_results, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Optional: Print a message indicating completion
cat("Merging complete. Results saved to", output_path, "\n")
