# Define the data
benign_control_variants_in_EM_genes <- c(1411, 756, 67, 168, 793, 477, 45, 297, 89, 169, 132, 222, 203, 543, 29, 213)
benign_control_variants_in_non_EM_genes <- c(24712, 10221, 1574, 2782, 8483, 6451, 797, 4447, 1961, 2414, 2300, 2242, 2688, 9168, 721, 5003)
pathogenic_variants_in_EM_genes <- c(349, 459, 123, 133, 261, 189, 432, 467, 101, 92, 42, 93, 394, 143, 35, 49)
pathogenic_variants_in_non_EM_genes <- c(7427, 3926, 1607, 2944, 2765, 2750, 433, 1796, 1058, 977, 289, 1115, 1791, 2405, 394, 1046)

categories <- c("Nonpolar_to_Nonpolar", "Nonpolar_to_Polar", "Nonpolar_to_Negatively_charged", "Nonpolar_to_Positively_charged",
                "Polar_to_Nonpolar", "Polar_to_Polar", "Polar_to_Negatively_charged", "Polar_to_Positively_charged",
                "Negatively_charged_to_Nonpolar", "Negatively_charged_to_Polar", "Negatively_charged_to_Negatively_charged",
                "Negatively_charged_to_Positively_charged", "Positively_charged_to_Nonpolar", "Positively_charged_to_Polar",
                "Positively_charged_to_Negatively_charged", "Positively_charged_to_Positively_charged")

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

# Manually create each contingency table and perform Fisher's exact test
cat("Fisher's Exact Test Results:\n")

for (i in 1:length(categories)) {
  cat("\nCategory: ", categories[i], "\n", sep="")
  
  # Create contingency table
  contingency_table <- data.frame(
    "pathogenic_variants" = c(pathogenic_variants_in_EM_genes[i], pathogenic_variants_in_non_EM_genes[i]),
    "control_variants" = c(benign_control_variants_in_EM_genes[i], benign_control_variants_in_non_EM_genes[i]),
    row.names = c("EM_genes", "non_EM_genes"),
    stringsAsFactors = FALSE
  )
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  cat("Contingency Table:\n")
  print(contingency_table)
  
  cat("P-Value: ", fisher_test$p.value, "\n")
  cat("Odds Ratio: ", fisher_test$estimate, "\n")
  cat("95% Confidence Interval: ", fisher_test$conf.int, "\n")
}

sink()  # Close the connection to the file

# Inform user that the results have been written
cat("Results have been written to", output_file, "\n")
