benign-epigenetic <- "
Nonpolar->Nonpolar 1411
Nonpolar->Polar 756
Nonpolar->Negatively_charged 67
Nonpolar->Positively_charged 168
Polar->Nonpolar 793
Polar->Polar 477
Polar->Negatively_charged 45
Polar->Positively_charged 297
Negatively_charged->Nonpolar 89
Negatively_charged->Polar 169
Negatively_charged->Negatively_charged 132
Negatively_charged->Positively_charged 222
Positively_charged->Nonpolar 203
Positively_charged->Polar 543
Positively_charged->Negatively_charged 29
Positively_charged->Positively_charged 213
"

benign-control <- "
Nonpolar->Nonpolar 24712
Nonpolar->Polar 10221
Nonpolar->Negatively_charged 1574
Nonpolar->Positively_charged 2782
Polar->Nonpolar 8483
Polar->Polar 6451
Polar->Negatively_charged 797
Polar->Positively_charged 4447
Negatively_charged->Nonpolar 1961
Negatively_charged->Polar 2414
Negatively_charged->Negatively_charged 2300
Negatively_charged->Positively_charged 2242
Positively_charged->Nonpolar 2688
Positively_charged->Polar 9168
Positively_charged->Negatively_charged 721
Positively_charged->Positively_charged 5003
"

pathogenic-epigenetic <- "
Nonpolar->Nonpolar 349
Nonpolar->Polar 459
Nonpolar->Negatively_charged 123
Nonpolar->Positively_charged 133
Polar->Nonpolar 261
Polar->Polar 189
Polar->Negatively_charged 432
Polar->Positively_charged 467
Negatively_charged->Nonpolar 101
Negatively_charged->Polar 92
Negatively_charged->Negatively_charged 42
Negatively_charged->Positively_charged 93
Positively_charged->Nonpolar 394
Positively_charged->Polar 143
Positively_charged->Negatively_charged 35
Positively_charged->Positively_charged 49
"

pathogenic-control <- "
Nonpolar->Nonpolar 7427
Nonpolar->Polar 3926
Nonpolar->Negatively_charged 1607
Nonpolar->Positively_charged 2944
Polar->Nonpolar 2765
Polar->Polar 2750
Polar->Negatively_charged 433
Polar->Positively_charged 1796
Negatively_charged->Nonpolar 1058
Negatively_charged->Polar 977
Negatively_charged->Negatively_charged 289
Negatively_charged->Positively_charged 1115
Positively_charged->Nonpolar 1791
Positively_charged->Polar 2405
Positively_charged->Negatively_charged 394
Positively_charged->Positively_charged 1046
"

# Define the contingency tables for each categorical change
benign_control_epigenetic <- c(1411, 756, 67, 168, 793, 477, 45, 297, 89, 169, 132, 222, 203, 543, 29, 213)
benign_control_non_epigenetic <- c(24712, 10221, 1574, 2782, 8483, 6451, 797, 4447, 1961, 2414, 2300, 2242, 2688, 9168, 721, 5003)
pathogenic_epigenetic <- c(349, 459, 123, 133, 261, 189, 432, 467, 101, 92, 42, 93, 394, 143, 35, 49)
pathogenic_non_epigenetic <- c(7427, 3926, 1607, 2944, 2765, 2750, 433, 1796, 1058, 977, 289, 1115, 1791, 2405, 394, 1046)

categories <- c("Nonpolar->Nonpolar", "Nonpolar->Polar", "Nonpolar->Negatively_charged", "Nonpolar->Positively_charged",
                "Polar->Nonpolar", "Polar->Polar", "Polar->Negatively_charged", "Polar->Positively_charged",
                "Negatively_charged->Nonpolar", "Negatively_charged->Polar", "Negatively_charged->Negatively_charged",
                "Negatively_charged->Positively_charged", "Positively_charged->Nonpolar", "Positively_charged->Polar",
                "Positively_charged->Negatively_charged", "Positively_charged->Positively_charged")

# Initialize an empty data frame to store results
results <- data.frame(Categorical_Change = character(), Pathogenic_Epigenetic = numeric(), Benign_Control_Epigenetic = numeric(),
                      Pathogenic_Non_Epigenetic = numeric(), Benign_Control_Non_Epigenetic = numeric(), P_Value = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through each category and perform Fisher's exact test
for (i in 1:length(categories)) {
  table <- matrix(c(pathogenic_epigenetic[i], benign_control_epigenetic[i],
                    pathogenic_non_epigenetic[i], benign_control_non_epigenetic[i]), nrow = 2)
  fisher_test <- fisher.test(table)
  p_value <- fisher_test$p.value
  results <- rbind(results, data.frame(Categorical_Change = categories[i],
                                       Pathogenic_Epigenetic = pathogenic_epigenetic[i],
                                       Benign_Control_Epigenetic = benign_control_epigenetic[i],
                                       Pathogenic_Non_Epigenetic = pathogenic_non_epigenetic[i],
                                       Benign_Control_Non_Epigenetic = benign_control_non_epigenetic[i],
                                       P_Value = p_value))
}

# Adjust for multiple comparisons using the Bonferroni method
results$Corrected_P_Value <- p.adjust(results$P_Value, method = "bonferroni")

# Write the results to a file
output_file <- "fisher_test_results.txt"
sink(output_file)

# Print original datasets
cat("Original Datasets:\n")
cat("Benign Control Epigenetic:\n")
print(benign_control_epigenetic)
cat("\nBenign Control Non-Epigenetic:\n")
print(benign_control_non_epigenetic)
cat("\nPathogenic Epigenetic:\n")
print(pathogenic_epigenetic)
cat("\nPathogenic Non-Epigenetic:\n")
print(pathogenic_non_epigenetic)
cat("\n\n")

# Print results
cat("Fisher's Exact Test Results:\n")
print(results)

sink()  # Close the connection to the file

# Inform user that the results have been written
cat("Results have been written to", output_file, "\n")

"""
contingency_table <- data.frame(
  "EM genes" = c(432, 45),
  "non-EM genes" = c(433, 797),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )
"""
