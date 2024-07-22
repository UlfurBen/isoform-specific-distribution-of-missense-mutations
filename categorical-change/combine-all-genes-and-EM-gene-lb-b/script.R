# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the data
dataset1 <- "
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

dataset2 <- "
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
