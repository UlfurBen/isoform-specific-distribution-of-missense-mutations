# Load necessary libraries
library(ggplot2)

# Read the data from enrichment.txt file
enrichment_data <- read.csv("enrichment.txt", header = TRUE)

# Identify main isoform (without hyphen in identifier)
enrichment_data$Is_Rare_Isoform <- grepl("-", enrichment_data$Identifier)

# Calculate the total number of variants in each category using only 'Variation count'
total_all_variants <- sum(enrichment_data$Variation.count)
total_main_variants <- sum(enrichment_data$Variation.count[!enrichment_data$Is_Rare_Isoform])
total_rare_variants <- sum(enrichment_data$Variation.count[enrichment_data$Is_Rare_Isoform])

# Calculate the percentages
percent_main_isoform <- (total_main_variants / total_all_variants) * 100
percent_all_isoforms <- 100  # Since this is the total of all variants
percent_rare_isoform <- (total_rare_variants / total_all_variants) * 100

# Prepare data for plotting
plot_data <- data.frame(
  Isoform_Category = c("Main Isoform", "All Isoforms", "Rare Isoforms"),
  Percentage = c(percent_main_isoform, percent_all_isoforms, percent_rare_isoform)
)

# Plot the data
ggplot(plot_data, aes(x = Isoform_Category, y = Percentage, fill = Isoform_Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Variation Distribution Across Isoforms", x = "Isoform Category", y = "Percentage of Variations (%)") +
  theme_minimal()
