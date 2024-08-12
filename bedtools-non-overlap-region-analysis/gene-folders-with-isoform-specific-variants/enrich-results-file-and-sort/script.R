# Set the file paths
input_file <- "~/results/filtered_results.bed"
output_dir <- "~/final-results"
output_file <- file.path(output_dir, "results.bed")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the BED file into a dataframe
bed_data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)

# Calculate the enrichment value and add it as the 7th column
bed_data$V7 <- bed_data$V6 / (bed_data$V3 - bed_data$V2)

# Sort the dataframe by the enrichment value (7th column)
sorted_bed_data <- bed_data[order(-bed_data$V7),]

# Save the sorted dataframe to a new BED file
write.table(sorted_bed_data, output_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

cat("File has been saved to:", output_file, "\n")
