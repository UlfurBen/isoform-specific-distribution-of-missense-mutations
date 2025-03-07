# Initialize an empty list to store results
results_list <- list()

# Define a function to perform Fisher's test and save results
perform_fisher_test <- function(data, change_name) {
  fisher_test <- fisher.test(data)
  result <- data.frame(
    Categorical_Change = change_name,
    Pathogenic_Epigenetic = data["pathogenic_variants", "EM_genes"],
    Benign_Control_Epigenetic = data["control_variants", "EM_genes"],
    Pathogenic_Non_Epigenetic = data["pathogenic_variants", "non-EM_genes"],
    Benign_Control_Non_Epigenetic = data["control_variants", "non-EM_genes"],
    P_Value = fisher_test$p.value,
    Odds_Ratio = fisher_test$estimate
  )
  return(result)
}

# List of contingency tables and their names
contingency_tables <- list(
  Nonpolar_to_Nonpolar = data.frame(
    "EM_genes" = c(349, 1411),
    "non-EM_genes" = c(7427, 24712),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Nonpolar_to_Polar = data.frame(
    "EM_genes" = c(459, 756),
    "non-EM_genes" = c(3926, 10221),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Nonpolar_to_Negatively_charged = data.frame(
    "EM_genes" = c(123, 67),
    "non-EM_genes" = c(1607, 1574),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Nonpolar_to_Positively_charged = data.frame(
    "EM_genes" = c(133, 168),
    "non-EM_genes" = c(2944, 2782),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Polar_to_Nonpolar = data.frame(
    "EM_genes" = c(261, 793),
    "non-EM_genes" = c(2765, 8483),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Polar_to_Polar = data.frame(
    "EM_genes" = c(189, 477),
    "non-EM_genes" = c(2750, 6451),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Polar_to_Negatively_charged = data.frame(
    "EM_genes" = c(432, 45),
    "non-EM_genes" = c(433, 797),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Polar_to_Positively_charged = data.frame(
    "EM_genes" = c(467, 297),
    "non-EM_genes" = c(1796, 4447),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Negatively_charged_to_Nonpolar = data.frame(
    "EM_genes" = c(101, 89),
    "non-EM_genes" = c(1058, 1961),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Negatively_charged_to_Polar = data.frame(
    "EM_genes" = c(92, 169),
    "non-EM_genes" = c(977, 2414),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Negatively_charged_to_Negatively_charged = data.frame(
    "EM_genes" = c(42, 132),
    "non-EM_genes" = c(289, 2300),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Negatively_charged_to_Positively_charged = data.frame(
    "EM_genes" = c(93, 222),
    "non-EM_genes" = c(1115, 2242),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Positively_charged_to_Nonpolar = data.frame(
    "EM_genes" = c(394, 203),
    "non-EM_genes" = c(1791, 2688),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Positively_charged_to_Polar = data.frame(
    "EM_genes" = c(143, 543),
    "non-EM_genes" = c(2405, 9168),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Positively_charged_to_Negatively_charged = data.frame(
    "EM_genes" = c(35, 29),
    "non-EM_genes" = c(394, 721),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  ),
  Positively_charged_to_Positively_charged = data.frame(
    "EM_genes" = c(49, 213),
    "non-EM_genes" = c(1046, 5003),
    row.names = c("pathogenic_variants", "control_variants"),
    stringsAsFactors = FALSE
  )
)

# Loop through each contingency table and perform Fisher's test
for (change_name in names(contingency_tables)) {
  results_list[[change_name]] <- perform_fisher_test(contingency_tables[[change_name]], change_name)
}

# Combine all results into one data frame
results <- do.call(rbind, results_list)

# Adjust for multiple comparisons using the Bonferroni method
results$Corrected_P_Value <- p.adjust(results$P_Value, method = "bonferroni")

# Write the results to a file
output_file <- "fisher_test_results.txt"
sink(output_file)

# Print original datasets
cat("Original Datasets:\n")
cat("Benign Control Epigenetic:\n")
print(benign_control_variants_in_EM_genes)
cat("\nBenign Control Non-Epigenetic:\n")
print(benign_control_variants_in_non_EM_genes)
cat("\nPathogenic Epigenetic:\n")
print(pathogenic_variants_in_EM_genes)
cat("\nPathogenic Non-Epigenetic:\n")
print(pathogenic_variants_in_non_EM_genes)
cat("\n\n")

# Print results
cat("Fisher's Exact Test Results:\n")
print(results)

sink()  # Close the connection to the file

# Inform user that the results have been written
cat("Results have been written to", output_file, "\n")
