# Define the output file
output_file <- "fisher-result.txt"

# Create all contingency tables
Nonpolar_to_Nonpolar <- data.frame(
  "EM_genes" = c(349, 1411),
  "non-EM_genes" = c(7427, 24712),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Nonpolar_to_Polar <- data.frame(
  "EM_genes" = c(459, 756),
  "non-EM_genes" = c(3926, 10221),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Nonpolar_to_Negatively_charged <- data.frame(
  "EM_genes" = c(123, 67),
  "non-EM_genes" = c(1607, 1574),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Nonpolar_to_Positively_charged <- data.frame(
  "EM_genes" = c(133, 168),
  "non-EM_genes" = c(2944, 2782),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Polar_to_Nonpolar <- data.frame(
  "EM_genes" = c(261, 793),
  "non-EM_genes" = c(2765, 8483),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Polar_to_Polar <- data.frame(
  "EM_genes" = c(189, 477),
  "non-EM_genes" = c(2750, 6451),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Polar_to_Negatively_charged <- data.frame(
  "EM_genes" = c(432, 45),
  "non-EM_genes" = c(433, 797),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Polar_to_Positively_charged <- data.frame(
  "EM_genes" = c(467, 297),
  "non-EM_genes" = c(1796, 4447),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Negatively_charged_to_Nonpolar <- data.frame(
  "EM_genes" = c(101, 89),
  "non-EM_genes" = c(1058, 1961),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Negatively_charged_to_Polar <- data.frame(
  "EM_genes" = c(92, 169),
  "non-EM_genes" = c(977, 2414),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Negatively_charged_to_Negatively_charged <- data.frame(
  "EM_genes" = c(42, 132),
  "non-EM_genes" = c(289, 2300),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Negatively_charged_to_Positively_charged <- data.frame(
  "EM_genes" = c(93, 222),
  "non-EM_genes" = c(1115, 2242),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Positively_charged_to_Nonpolar <- data.frame(
  "EM_genes" = c(394, 203),
  "non-EM_genes" = c(1791, 2688),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Positively_charged_to_Polar <- data.frame(
  "EM_genes" = c(143, 543),
  "non-EM_genes" = c(2405, 9168),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Positively_charged_to_Negatively_charged <- data.frame(
  "EM_genes" = c(35, 29),
  "non-EM_genes" = c(394, 721),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

Positively_charged_to_Positively_charged <- data.frame(
  "EM_genes" = c(49, 213),
  "non-EM_genes" = c(1046, 5003),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)

# Perform Fisher's exact test and write the results to the file
write_fisher_test_result <- function(table, description) {
  result <- fisher.test(table)
  write(paste(description, ":\n", capture.output(result), "\n"), file = output_file, append = TRUE)
}

# Clear the file before writing new results
cat("", file = output_file)

# List of tables and their descriptions
tables <- list(
  list(table = Nonpolar_to_Nonpolar, description = "Fisher test result for Nonpolar_to_Nonpolar"),
  list(table = Nonpolar_to_Polar, description = "Fisher test result for Nonpolar_to_Polar"),
  list(table = Nonpolar_to_Negatively_charged, description = "Fisher test result for Nonpolar_to_Negatively_charged"),
  list(table = Nonpolar_to_Positively_charged, description = "Fisher test result for Nonpolar_to_Positively_charged"),
  list(table = Polar_to_Nonpolar, description = "Fisher test result for Polar_to_Nonpolar"),
  list(table = Polar_to_Polar, description = "Fisher test result for Polar_to_Polar"),
  list(table = Polar_to_Negatively_charged, description = "Fisher test result for Polar_to_Negatively_charged"),
  list(table = Polar_to_Positively_charged, description = "Fisher test result for Polar_to_Positively_charged"),
  list(table = Negatively_charged_to_Nonpolar, description = "Fisher test result for Negatively_charged_to_Nonpolar"),
  list(table = Negatively_charged_to_Polar, description = "Fisher test result for Negatively_charged_to_Polar"),
  list(table = Negatively_charged_to_Negatively_charged, description = "Fisher test result for Negatively_charged_to_Negatively_charged"),
  list(table = Negatively_charged_to_Positively_charged, description = "Fisher test result for Negatively_charged_to_Positively_charged"),
  list(table = Positively_charged_to_Nonpolar, description = "Fisher test result for Positively_charged_to_Nonpolar"),
  list(table = Positively_charged_to_Polar, description = "Fisher test result for Positively_charged_to_Polar"),
  list(table = Positively_charged_to_Negatively_charged, description = "Fisher test result for Positively_charged_to_Negatively_charged"),
  list(table = Positively_charged_to_Positively_charged, description = "Fisher test result for Positively_charged_to_Positively_charged")
)

# Perform tests and write results
for (item in tables) {
  write_fisher_test_result(item$table, item$description)
}

# Inform user that the results have been written
cat("Results have been written to", output_file, "\n")
