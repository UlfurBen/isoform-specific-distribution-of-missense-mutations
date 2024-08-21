# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

# Function to read and process data from a file
process_file_data <- function(file_path, label) {
  # Read the file
  df <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Set column names
  colnames(df) <- c("Category", "Count")
  # Add a label to identify the dataset
  df$Dataset <- label
  return(df)
}

# File paths and labels
file1 <- "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_EM_genes_aa_change_count.txt"
file2 <- "homo_sapiens_variation_missense_ClinVar_pathogenic_likely_pathogenic_only_aa_change_count.txt"

label1 <- "EM genes"
label2 <- "All genes"

# Process the files
df1 <- process_file_data(file1, label1)
df2 <- process_file_data(file2, label2)

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

# Create the stacked bar plot with ggplot2
p <- ggplot(combined_df, aes(x = Category, y = Percentage, fill = Dataset, 
                             text = paste("Category:", Category, "<br>Count:", Count, "<br>Percentage:", round(Percentage, 2), "%"))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage of Each Categorical Change of ClinVar Missense Variants (Pathogenic/Likely Pathogenic)", 
       y = "Percentage", 
       x = "Categorical Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

# Convert ggplot to an interactive plotly plot
interactive_plot <- ggplotly(p, tooltip = "text")

# Show the interactive plot
interactive_plot
