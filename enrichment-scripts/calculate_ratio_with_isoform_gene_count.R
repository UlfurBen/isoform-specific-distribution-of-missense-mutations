input_file <- "appended_with_isoform_counts.txt"
output_file <- "calculate_ratio_with_isoform_gene_count.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Calculate the ratio of Count to Length
data <- transform(data, Ratio = as.numeric(Count) / (as.numeric(Isoform_Count)*as.numeric(Length)))

# Write the output file
write.table(data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print a message indicating completion
cat("Ratios have been calculated and saved to", output_file, "\n")
