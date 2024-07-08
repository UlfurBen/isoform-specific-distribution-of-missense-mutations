# Load required libraries
library(ggplot2)
library(reshape2)
library(ggiraph)

# Define file paths
input_file <- "enrichment.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Calculate ratio for each relevant pathogenicity type divided by variant_space
data$Likely_benign_Ratio <- data$Likely.benign / data$variant_space
data$Uncertain_Ratio <- data$uncertain / data$variant_space
data$Benign_Ratio <- data$Benign / data$variant_space
data$Pathogenic_Ratio <- data$Pathogenic / data$variant_space
data$Likely_pathogenic_Ratio <- data$Likely.pathogenic / data$variant_space

# Calculate total counts for each identifier
data$TotalCount <- rowSums(data[, c("Likely_benign_Ratio", "Uncertain_Ratio", "Benign_Ratio", "Pathogenic_Ratio", "Likely_pathogenic_Ratio")])

# Sort data by TotalCount in descending order
data <- data[order(-data$TotalCount), ]

# Melt the data for ggplot2 to plot ratios
melted_data <- melt(data, id.vars = c("Identifier", "TotalCount", "variant_space"), 
                    measure.vars = c("Likely_benign_Ratio", "Uncertain_Ratio", "Benign_Ratio", "Pathogenic_Ratio", "Likely_pathogenic_Ratio"), 
                    variable.name = "Category", value.name = "Ratio")

# Create the stacked bar plot with ggiraph
p <- ggplot(melted_data, aes(x = reorder(Identifier, -TotalCount), y = Ratio, fill = Category)) +
  geom_bar_interactive(stat = "identity", aes(tooltip = Identifier), width = 0.8) +
  labs(x = "Identifier", y = "Ratio (Variant Count / Sequence Length)", title = "Ratio of Variant Counts by Category for Each Identifier") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Display the plot
girafe(ggobj = p, width_svg = 10)

# Print completion message
cat("Plot displayed.\n")
