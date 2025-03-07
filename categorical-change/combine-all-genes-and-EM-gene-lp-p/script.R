# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the data
dataset1 <- "
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

dataset2 <- "
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

# Convert the datasets into data frames
process_data <- function(data, label) {
  df <- read.table(text = data, header = FALSE, sep = " ")
  colnames(df) <- c("Category", "Count")
  df$Dataset <- label
  return(df)
}

df1 <- process_data(dataset1, "EM genes")
df2 <- process_data(dataset2, "all genes")

# Combine the datasets
combined_df <- bind_rows(df1, df2)

# Calculate the percentages
combined_df <- combined_df %>%
  group_by(Dataset) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Calculate the total percentage for ordering
total_percentage_df <- combined_df %>%
  group_by(Category) %>%
  summarise(TotalPercentage = sum(Percentage))

# Merge total percentages back to the combined data
combined_df <- combined_df %>%
  left_join(total_percentage_df, by = "Category")

# Order the data by total percentage in descending order
combined_df <- combined_df %>%
  mutate(Category = factor(Category, levels = total_percentage_df$Category[order(-total_percentage_df$TotalPercentage)]))

# Create the stacked bar plot
ggplot(combined_df, aes(x = Category, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage of Each Categorical Change of ClinVar Missense Variants", 
       y = "Percentage", 
       x = "Categorical Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")
