# Define file paths
input_file <- "enrichment.txt"
fasta_file <- "uniprot_sprot_varsplic.fasta"
output_file <- "enrichment_filtered_with_seq_len_more_than_before.txt"

# Read the input file
input_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Function to read FASTA file
read_fasta <- function(file) {
  lines <- readLines(file)
  fasta_sequences <- list()
  seq_name <- NULL
  seq <- NULL
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(seq_name)) {
        fasta_sequences[[seq_name]] <- paste(seq, collapse = "")
      }
      seq_name <- sub(" .*", "", sub("^>", "", line))
      seq <- c()
    } else {
      seq <- c(seq, line)
    }
  }
  if (!is.null(seq_name)) {
    fasta_sequences[[seq_name]] <- paste(seq, collapse = "")
  }
  return(fasta_sequences)
}

# Function to get sequence length from UniProt REST API
get_sequence_length_from_uniprot <- function(isoform_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", isoform_id, ".fasta")
  response <- tryCatch(
    {
      readLines(url)
    },
    error = function(e) {
      message(paste("Error fetching data for", isoform_id))
      return(NULL)
    }
  )
  if (!is.null(response)) {
    seq <- paste(response[!startsWith(response, ">")], collapse = "")
    return(nchar(seq))
  } else {
    return(NA)
  }
}

# Read the FASTA file
fasta_sequences <- read_fasta(fasta_file)

# Initialize lists to store lengths for isoforms
isoform_lengths <- list()

# Process each isoform in the input data
for (i in 1:nrow(input_data)) {
  isoform_id <- input_data[i, 2]  # Assuming the second column contains the identifier
  length <- NA
  
  # Check if the isoform is in the FASTA file and if the line contains '_HUMAN'
  if (isoform_id %in% names(fasta_sequences)) {
    if (grepl("_HUMAN", fasta_sequences[[isoform_id]])) {
      length <- nchar(fasta_sequences[[isoform_id]])
      print(paste("Found in FASTA and contains _HUMAN:", isoform_id, "Length:", length))  # Debug print statement
    } else {
      print(paste("Found in FASTA but does not contain _HUMAN:", isoform_id))  # Debug print statement
    }
  } else {
    # Fetch from UniProt REST API
    length <- get_sequence_length_from_uniprot(isoform_id)
    if (!is.na(length)) {
      print(paste("Fetched from UniProt API:", isoform_id, "Length:", length))  # Debug print statement
    } else {
      print(paste("Not found in FASTA and not fetched from UniProt API:", isoform_id))  # Debug print statement
    }
  }
  
  # Store the length in the list
  isoform_lengths[[isoform_id]] <- length
}

# Add the lengths to the input data
input_data$seq_len <- sapply(input_data[[2]], function(id) isoform_lengths[[id]])

# Write the merged results to the final output file
write.table(input_data, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print message indicating completion
cat("Isoform lengths have been processed and appended to", output_file, "\n")
