# Load the required library
library(ggplot2)

# Define file paths
input_file <- "merged_calculate_missense_variant_enrichment_ratio_with_gene_total_variant_count.txt"

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Extract the gene identifier (part before the hyphen) from the Isoform column
data$Gene <- sub("-.*", "", data$Isoform)

# Get unique genes
unique_genes <- unique(data$Gene)

# Create a plot for each gene and save as PDF if it contains more than one data point
for (gene in unique_genes) {
  # Subset data for the current gene
  gene_data <- subset(data, Gene == gene)
  
  # Check if there is more than one data point
  if (nrow(gene_data) > 1) {
    # Determine which isoform is the main one (without hyphen)
    gene_data$IsMainIsoform <- ifelse(grepl("-", gene_data$Isoform), FALSE, TRUE)
    
    # Create the plot
    p <- ggplot(gene_data, aes(x = Isoform, y = Ratio)) +
      geom_point(aes(color = IsMainIsoform), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
      labs(title = paste("Dot Plot for Gene:", gene), x = "Isoform", y = "Ratio") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save the plot to a PDF file
    pdf(file = paste0(gene, "_dot_plot.pdf"))
    print(p)
    dev.off()
  }
}
