# Go to url https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/ and download
# homo_sapiens_variation.txt.gz and gunzip

# run command grep -i 'missense' homo_sapiens_variation.txt | grep -i 'ClinVar' | grep -vi 'somatic' > homo_sapiens_variation_missense_ClinVar_not_somatic.txt


# Rename columns as necessary

# Load necessary libraries
library(dplyr)
library(tidyr)

# Define the input and output file paths
input_file <- "homo_sapiens_variation_missense_ClinVar_not_somatic.txt"
output_file <- "homo_sapiens_variation_missense_ClinVar_not_somatic.bed"

# Define the column names manually
column_names <- c("GeneName", "AC", "VariantAAChange", "SourceDBID", "ConsequenceType", "ClinicalSignificance",
                  "PhenotypeDisease", "PhenotypeDiseaseSource", "CytogeneticBand", "ChromosomeCoordinate", 
                  "EnsemblGeneID", "EnsemblTranscriptID", "EnsemblTranslationID", "Evidence")

# Read the input file without headers and assign the column names
data <- read.delim(input_file, header = FALSE, sep = "\t", col.names = column_names)

# Select and rename columns as per the specification
output_data <- data %>%
  select(GeneName = GeneName,
         AC = AC,
         VariantAAChange = VariantAAChange,
         SourceDBID = SourceDBID,
         ConsequenceType = ConsequenceType,
         ClinicalSignificance = ClinicalSignificance,
         PhenotypeDisease = PhenotypeDisease,
         PhenotypeDiseaseSource = PhenotypeDiseaseSource,
         CytogeneticBand = CytogeneticBand,
         ChromosomeCoordinate = ChromosomeCoordinate,
         EnsemblGeneID = EnsemblGeneID,
         EnsemblTranscriptID = EnsemblTranscriptID,
         EnsemblTranslationID = EnsemblTranslationID,
         Evidence = Evidence)

# Write the output to a BED file
write.table(output_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("BED file has been written to", output_file, "\n")







# Get chr, start and end values

# Load necessary libraries
library(dplyr)
library(stringr)

# Define the input and output file paths
input_file <- "homo_sapiens_variation_missense_ClinVar_not_somatic.bed"
output_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"

# Read the input file
data <- read.delim(input_file, header = TRUE, sep = "\t")

# Create new columns chr, chromStart, chromEnd, and aachange
data <- data %>%
  mutate(
    chr = str_extract(ChromosomeCoordinate, "(?<=NC_)(0*)([1-9][0-9]*)") %>% str_extract("[1-9][0-9]*"),
    chromStart = as.numeric(str_extract(ChromosomeCoordinate, "(?<=:g\\.)([0-9]+)")),
    chromEnd = chromStart,
    aachange = str_extract(ChromosomeCoordinate, "(?<=[0-9])[A-Z]>[A-Z]")
  ) %>%
  # Reorder columns: place chr, chromStart, chromEnd, aachange in front of all original columns
  select(chr, chromStart, chromEnd, aachange, everything())

# Write the output to a new BED file
write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Filtered and modified data have been written to", output_file, "\n")








# Final filtering steps before use in enrichment analysis

#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(dplyr)

# File paths
input_file <- "homo_sapiens_variation_missense_ClinVar_filtered_relevancy.bed"
output_file <- "variants_benign_pathogenic_non_vus_non_conflicting_including_non_annotation.bed"

# Read the BED file and skip the header row
bed_data <- fread(input_file, skip = 1, header = FALSE)

# Define the exclusion terms
exclusion_terms <- c("variant of uncertain significance", "conflicting")

# Exclude rows with conflicting terms in the same column
filtered_data <- bed_data %>%
  filter(!grepl(paste(exclusion_terms, collapse = "|"), V10, ignore.case = TRUE))

# Write the filtered data to a new file
fwrite(filtered_data, output_file, sep = "\t", col.names = FALSE)
cat("Filtered file has been processed and saved as", output_file, "\n")
