dplyr_path <- "/hpchome/ubf2/R/x86_64-pc-linux-gnu-library/4.1/dplyr"
tidyr_path <- "/hpcapps/lib-mimir/software/R/4.1.2-foss-2021b/lib64/R/library/tidyr"
stringr_path <- "/hpcapps/lib-mimir/software/R/4.1.2-foss-2021b/lib64/R/library/stringr"

.libPaths(c(dirname(dplyr_path), .libPaths()))
.libPaths(c(dirname(tidyr_path), .libPaths()))
.libPaths(c(dirname(stringr_path), .libPaths()))

# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)

# Define output files
output_bed_file <- "Homo_sapiens.GRCh37.87.bed"
output_bed_file_with_headers <- "Homo_sapiens.GRCh37.87_with_headers.bed"

# Read the GTF file, skipping lines starting with '#!'
gtf_file <- "Homo_sapiens.GRCh37.87.gtf"
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", quote="")

# Assign column names
colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Function to extract attribute values into a named list
extract_attributes <- function(attribute_string) {
  attribute_list <- str_split(attribute_string, "; ")[[1]]
  attribute_list <- attribute_list[attribute_list != ""]
  attribute_pairs <- str_split(attribute_list, " ")
  attribute_names <- sapply(attribute_pairs, function(x) x[1])
  attribute_values <- sapply(attribute_pairs, function(x) gsub('"', '', x[2]))
  names(attribute_values) <- attribute_names
  return(as.list(attribute_values))
}

# Apply the function to the attribute column
attributes_list <- lapply(gtf_data$attribute, extract_attributes)

# Combine the list of attributes into a data frame
attributes_df <- do.call(bind_rows, attributes_list)

# Ensure all columns are present in attributes_df
expected_columns <- c("gene_id", "gene_version", "transcript_id", "transcript_version", "exon_number", 
                      "gene_name", "gene_source", "gene_biotype", "transcript_name", "transcript_source", 
                      "transcript_biotype", "protein_id", "protein_version", "tag")
for (col in expected_columns) {
  if (!(col %in% colnames(attributes_df))) {
    attributes_df[[col]] <- NA
  }
}

# Combine the attributes data frame with the original GTF data frame
bed_df <- cbind(gtf_data %>% select(-attribute), attributes_df)

# Write to a BED file without headers
write.table(bed_df, file=output_bed_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Read the output BED file
bed_df <- read.table(output_bed_file, sep="\t", header=FALSE)

# Set column names
colnames(bed_df) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", expected_columns)

# Display the first 10 rows to confirm changes
print(head(bed_df, 10))

# Write the modified data frame back to a new file with headers
write.table(bed_df, file=output_bed_file_with_headers, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
