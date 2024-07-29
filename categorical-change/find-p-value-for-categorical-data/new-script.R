# Define the contingency tables for each categorical change
benign_control_variants_in_EM_genes <- c(1411, 756, 67, 168, 793, 477, 45, 297, 89, 169, 132, 222, 203, 543, 29, 213)
benign_control_variants_in_non_EM_genes <- c(24712, 10221, 1574, 2782, 8483, 6451, 797, 4447, 1961, 2414, 2300, 2242, 2688, 9168, 721, 5003)
pathogenic_variants_in_EM_genes <- c(349, 459, 123, 133, 261, 189, 432, 467, 101, 92, 42, 93, 394, 143, 35, 49)
pathogenic_variants_in_non_EM_genes <- c(7427, 3926, 1607, 2944, 2765, 2750, 433, 1796, 1058, 977, 289, 1115, 1791, 2405, 394, 1046)

categories <- c("Nonpolar->Nonpolar", "Nonpolar->Polar", "Nonpolar->Negatively_charged", "Nonpolar->Positively_charged",
                "Polar->Nonpolar", "Polar->Polar", "Polar->Negatively_charged", "Polar->Positively_charged",
                "Negatively_charged->Nonpolar", "Negatively_charged->Polar", "Negatively_charged->Negatively_charged",
                "Negatively_charged->Positively_charged", "Positively_charged->Nonpolar", "Positively_charged->Polar",
                "Positively_charged->Negatively_charged", "Positively_charged->Positively_charged")

# Initialize an empty data frame to store results
results <- data.frame(Categorical_Change = character(), Pathogenic_Epigenetic = numeric(), Benign_Control_Epigenetic = numeric(),
                      Pathogenic_Non_Epigenetic = numeric(), Benign_Control_Non_Epigenetic = numeric(), P_Value = numeric(), 
                      Odds_Ratio = numeric(), stringsAsFactors = FALSE)

# Loop through each category and perform Fisher's exact test
for (i in 1:length(categories)) {
  table <- matrix(c(pathogenic_variants_in_EM_genes[i], benign_control_variants_in_EM_genes[i],
                    pathogenic_variants_in_non_EM_genes[i], benign_control_variants_in_the__non_EM_genes[i]), nrow = 2)
  fisher_test <- fisher.test(table)
  p_value <- fisher_test$p.value
  odds_ratio <- fisher_test$estimate
  results <- rbind(results, data.frame(Categorical_Change = categories[i],
                                       Pathogenic_Epigenetic = pathogenic_variants_in_EM_genes[i],
                                       Benign_Control_Epigenetic = benign_control_variants_in_EM_genes[i],
                                       Pathogenic_Non_Epigenetic = pathogenic_variants_in_non_EM_genes[i],
                                       Benign_Control_Non_Epigenetic = benign_control_variants_in_non_EM_genes[i],
                                       P_Value = p_value,
                                       Odds_Ratio = odds_ratio))
}

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
