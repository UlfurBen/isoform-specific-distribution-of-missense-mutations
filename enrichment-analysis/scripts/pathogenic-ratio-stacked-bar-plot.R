# Load required libraries
library(ggplot2)
library(ggiraph)

# Define file paths
input_file <- "enrichment.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Calculate the ratio of pathogenic mutation count divided by variant_space
data$Pathogenic_Ratio <- data$Pathogenic / data$variant_space

# Sort data by Pathogenic_Ratio in descending order
data <- data[order(-data$Pathogenic_Ratio), ]

# Create the bar plot with ggiraph
p <- ggplot(data, aes(x = reorder(Identifier, -Pathogenic_Ratio), y = Pathogenic_Ratio)) +
  geom_bar_interactive(stat = "identity", aes(tooltip = Identifier, fill = "Pathogenic_Ratio"), width = 0.8) +
  labs(x = "Identifier", y = "Pathogenic Mutation Count / Sequence Length", title = "Pathogenic Mutation Ratio for Each Identifier") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Display the plot
girafe(ggobj = p, width_svg = 10)

# Print completion message
cat("Plot displayed.\n")
