# Load required libraries
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(httr)
library(jsonlite)

# Load the data from gene_names.csv file
data <- read.csv("important_genes.csv", stringsAsFactors = FALSE, header=TRUE)

# Extract the Gene_Name column
gene_names <- data[[1]]  # Assuming the single column contains gene names

# Function to retrieve the canonical ENST value for a given gene symbol
get_canonical_enst <- function(gene_symbol) {
  url <- paste0("https://rest.ensembl.org/lookup/symbol/homo_sapiens/", gene_symbol)
  response <- GET(url, add_headers("Content-Type" = "application/json"))
  
  if (status_code(response) == 200) {
    gene_data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    if (!is.null(gene_data$canonical_transcript)) {
      # Extract the canonical_transcript and remove the version number after the period
      canonical_transcript <- strsplit(gene_data$canonical_transcript, "\\.")[[1]][1]
      return(canonical_transcript)
    }
  } else {
    message("Failed to retrieve data for gene symbol: ", gene_symbol)
  }
  return(NA)  # Return NA if no canonical transcript is found or if the request fails
}

# Apply the function to each gene name to get the canonical ENST values
enst_values <- sapply(gene_names, get_canonical_enst)

# Append the canonical ENST values to the data frame, only if they are not NA
data$Canonical_ENST <- enst_values

# Filter out rows where Canonical_ENST is NA
data <- data[!is.na(data$Canonical_ENST), ]

# Output the result to a new CSV file with the appended canonical ENST values
output_file <- "important_genes_with_enst.csv"
write.csv(data, file = output_file, row.names = FALSE, quote = FALSE)

cat("Canonical ENST values have been appended and saved to", output_file, "\n")
