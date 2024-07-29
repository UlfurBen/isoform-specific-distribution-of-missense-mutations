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

# Initialize an empty data frame to store results
results <- data.frame(Categorical_Change = character(), Pathogenic_Epigenetic = numeric(), Benign_Control_Epigenetic = numeric(),
                      Pathogenic_Non_Epigenetic = numeric(), Benign_Control_Non_Epigenetic = numeric(), P_Value = numeric(), 
                      Odds_Ratio = numeric(), stringsAsFactors = FALSE)

# Manually create each contingency table and perform Fisher's exact test

# Nonpolar_to_Nonpolar
Nonpolar_to_Nonpolar <- data.frame(
  "EM_genes" = c(349, 1411),
  "non-EM_genes" = c(7427, 24712),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Nonpolar_to_Nonpolar)
results <- rbind(results, data.frame(
  Categorical_Change = "Nonpolar_to_Nonpolar",
  Pathogenic_Epigenetic = 349,
  Benign_Control_Epigenetic = 1411,
  Pathogenic_Non_Epigenetic = 7427,
  Benign_Control_Non_Epigenetic = 24712,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Nonpolar_to_Polar
Nonpolar_to_Polar <- data.frame(
  "EM_genes" = c(459, 756),
  "non-EM_genes" = c(3926, 10221),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Nonpolar_to_Polar)
results <- rbind(results, data.frame(
  Categorical_Change = "Nonpolar_to_Polar",
  Pathogenic_Epigenetic = 459,
  Benign_Control_Epigenetic = 756,
  Pathogenic_Non_Epigenetic = 3926,
  Benign_Control_Non_Epigenetic = 10221,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Nonpolar_to_Negatively_charged
Nonpolar_to_Negatively_charged <- data.frame(
  "EM_genes" = c(123, 67),
  "non-EM_genes" = c(1607, 1574),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Nonpolar_to_Negatively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Nonpolar_to_Negatively_charged",
  Pathogenic_Epigenetic = 123,
  Benign_Control_Epigenetic = 67,
  Pathogenic_Non_Epigenetic = 1607,
  Benign_Control_Non_Epigenetic = 1574,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Nonpolar_to_Positively_charged
Nonpolar_to_Positively_charged <- data.frame(
  "EM_genes" = c(133, 168),
  "non-EM_genes" = c(2944, 2782),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Nonpolar_to_Positively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Nonpolar_to_Positively_charged",
  Pathogenic_Epigenetic = 133,
  Benign_Control_Epigenetic = 168,
  Pathogenic_Non_Epigenetic = 2944,
  Benign_Control_Non_Epigenetic = 2782,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Polar_to_Nonpolar
Polar_to_Nonpolar <- data.frame(
  "EM_genes" = c(261, 793),
  "non-EM_genes" = c(2765, 8483),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Polar_to_Nonpolar)
results <- rbind(results, data.frame(
  Categorical_Change = "Polar_to_Nonpolar",
  Pathogenic_Epigenetic = 261,
  Benign_Control_Epigenetic = 793,
  Pathogenic_Non_Epigenetic = 2765,
  Benign_Control_Non_Epigenetic = 8483,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Polar_to_Polar
Polar_to_Polar <- data.frame(
  "EM_genes" = c(189, 477),
  "non-EM_genes" = c(2750, 6451),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Polar_to_Polar)
results <- rbind(results, data.frame(
  Categorical_Change = "Polar_to_Polar",
  Pathogenic_Epigenetic = 189,
  Benign_Control_Epigenetic = 477,
  Pathogenic_Non_Epigenetic = 2750,
  Benign_Control_Non_Epigenetic = 6451,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Polar_to_Negatively_charged
Polar_to_Negatively_charged <- data.frame(
  "EM_genes" = c(432, 45),
  "non-EM_genes" = c(433, 797),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Polar_to_Negatively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Polar_to_Negatively_charged",
  Pathogenic_Epigenetic = 432,
  Benign_Control_Epigenetic = 45,
  Pathogenic_Non_Epigenetic = 433,
  Benign_Control_Non_Epigenetic = 797,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Polar_to_Positively_charged
Polar_to_Positively_charged <- data.frame(
  "EM_genes" = c(467, 297),
  "non-EM_genes" = c(1796, 4447),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Polar_to_Positively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Polar_to_Positively_charged",
  Pathogenic_Epigenetic = 467,
  Benign_Control_Epigenetic = 297,
  Pathogenic_Non_Epigenetic = 1796,
  Benign_Control_Non_Epigenetic = 4447,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Negatively_charged_to_Nonpolar
Negatively_charged_to_Nonpolar <- data.frame(
  "EM_genes" = c(101, 89),
  "non-EM_genes" = c(1058, 1961),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Negatively_charged_to_Nonpolar)
results <- rbind(results, data.frame(
  Categorical_Change = "Negatively_charged_to_Nonpolar",
  Pathogenic_Epigenetic = 101,
  Benign_Control_Epigenetic = 89,
  Pathogenic_Non_Epigenetic = 1058,
  Benign_Control_Non_Epigenetic = 1961,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Negatively_charged_to_Polar
Negatively_charged_to_Polar <- data.frame(
  "EM_genes" = c(92, 169),
  "non-EM_genes" = c(977, 2414),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Negatively_charged_to_Polar)
results <- rbind(results, data.frame(
  Categorical_Change = "Negatively_charged_to_Polar",
  Pathogenic_Epigenetic = 92,
  Benign_Control_Epigenetic = 169,
  Pathogenic_Non_Epigenetic = 977,
  Benign_Control_Non_Epigenetic = 2414,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Negatively_charged_to_Negatively_charged
Negatively_charged_to_Negatively_charged <- data.frame(
  "EM_genes" = c(42, 132),
  "non-EM_genes" = c(289, 2300),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Negatively_charged_to_Negatively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Negatively_charged_to_Negatively_charged",
  Pathogenic_Epigenetic = 42,
  Benign_Control_Epigenetic = 132,
  Pathogenic_Non_Epigenetic = 289,
  Benign_Control_Non_Epigenetic = 2300,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Negatively_charged_to_Positively_charged
Negatively_charged_to_Positively_charged <- data.frame(
  "EM_genes" = c(93, 222),
  "non-EM_genes" = c(1115, 2242),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Negatively_charged_to_Positively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Negatively_charged_to_Positively_charged",
  Pathogenic_Epigenetic = 93,
  Benign_Control_Epigenetic = 222,
  Pathogenic_Non_Epigenetic = 1115,
  Benign_Control_Non_Epigenetic = 2242,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Positively_charged_to_Nonpolar
Positively_charged_to_Nonpolar <- data.frame(
  "EM_genes" = c(394, 203),
  "non-EM_genes" = c(1791, 2688),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Positively_charged_to_Nonpolar)
results <- rbind(results, data.frame(
  Categorical_Change = "Positively_charged_to_Nonpolar",
  Pathogenic_Epigenetic = 394,
  Benign_Control_Epigenetic = 203,
  Pathogenic_Non_Epigenetic = 1791,
  Benign_Control_Non_Epigenetic = 2688,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Positively_charged_to_Polar
Positively_charged_to_Polar <- data.frame(
  "EM_genes" = c(143, 543),
  "non-EM_genes" = c(2405, 9168),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Positively_charged_to_Polar)
results <- rbind(results, data.frame(
  Categorical_Change = "Positively_charged_to_Polar",
  Pathogenic_Epigenetic = 143,
  Benign_Control_Epigenetic = 543,
  Pathogenic_Non_Epigenetic = 2405,
  Benign_Control_Non_Epigenetic = 9168,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Positively_charged_to_Negatively_charged
Positively_charged_to_Negatively_charged <- data.frame(
  "EM_genes" = c(35, 29),
  "non-EM_genes" = c(394, 721),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Positively_charged_to_Negatively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Positively_charged_to_Negatively_charged",
  Pathogenic_Epigenetic = 35,
  Benign_Control_Epigenetic = 29,
  Pathogenic_Non_Epigenetic = 394,
  Benign_Control_Non_Epigenetic = 721,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

# Positively_charged_to_Positively_charged
Positively_charged_to_Positively_charged <- data.frame(
  "EM_genes" = c(49, 213),
  "non-EM_genes" = c(1046, 5003),
  row.names = c("pathogenic_variants", "control_variants"),
  stringsAsFactors = FALSE
)
fisher_test <- fisher.test(Positively_charged_to_Positively_charged)
results <- rbind(results, data.frame(
  Categorical_Change = "Positively_charged_to_Positively_charged",
  Pathogenic_Epigenetic = 49,
  Benign_Control_Epigenetic = 213,
  Pathogenic_Non_Epigenetic = 1046,
  Benign_Control_Non_Epigenetic = 5003,
  P_Value = fisher_test$p.value,
  Odds_Ratio = fisher_test$estimate
))

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
