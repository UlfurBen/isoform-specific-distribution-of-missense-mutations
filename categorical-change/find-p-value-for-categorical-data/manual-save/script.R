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

# Perform Fisher's exact test and write the results to the file with Bonferroni correction
write_fisher_test_result <- function(table, description, num_comparisons) {
  result <- fisher.test(table)
  p_value <- result$p.value
  bonferroni_p_value <- min(p_value * num_comparisons, 1)  # Apply Bonferroni correction and cap at 1
  result_string <- paste(
    description, ":\n",
    "Original p-value = ", format(p_value, digits = 5), "\n",
    "Bonferroni corrected p-value = ", format(bonferroni_p_value, digits = 5), "\n",
    "Odds ratio = ", format(result$estimate, digits = 5), "\n",
    "95% CI = [", format(result$conf.int[1], digits = 5), ", ", format(result$conf.int[2], digits = 5), "]\n",
    sep = ""
  )
  write(result_string, file = output_file, append = TRUE)
}

# Clear the file before writing new results
cat("", file = output_file)

# List of tables and their descriptions
tables <- list(
  list(table = Nonpolar_to_Nonpolar, description = "Nonpolar_to_Nonpolar"),
  list(table = Nonpolar_to_Polar, description = "Nonpolar_to_Polar"),
  list(table = Nonpolar_to_Negatively_charged, description = "Nonpolar_to_Negatively_charged"),
  list(table = Nonpolar_to_Positively_charged, description = "Nonpolar_to_Positively_charged"),
  list(table = Polar_to_Nonpolar, description = "Polar_to_Nonpolar"),
  list(table = Polar_to_Polar, description = "Polar_to_Polar"),
  list(table = Polar_to_Negatively_charged, description = "Polar_to_Negatively_charged"),
  list(table = Polar_to_Positively_charged, description = "Polar_to_Positively_charged"),
  list(table = Negatively_charged_to_Nonpolar, description = "Negatively_charged_to_Nonpolar"),
  list(table = Negatively_charged_to_Polar, description = "Negatively_charged_to_Polar"),
  list(table = Negatively_charged_to_Negatively_charged, description = "Negatively_charged_to_Negatively_charged"),
  list(table = Negatively_charged_to_Positively_charged, description = "Negatively_charged_to_Positively_charged"),
  list(table = Positively_charged_to_Nonpolar, description = "Positively_charged_to_Nonpolar"),
  list(table = Positively_charged_to_Polar, description = "Positively_charged_to_Polar"),
  list(table = Positively_charged_to_Negatively_charged, description = "Positively_charged_to_Negatively_charged"),
  list(table = Positively_charged_to_Positively_charged, description = "Positively_charged_to_Positively_charged")
)

# Count the number of comparisons
num_comparisons <- length(tables)

# Perform tests and write results with Bonferroni correction
for (item in tables) {
  write_fisher_test_result(item$table, item$description, num_comparisons)
}

# Inform user that the results have been written
cat("Results have been written to", output_file, "\n")
