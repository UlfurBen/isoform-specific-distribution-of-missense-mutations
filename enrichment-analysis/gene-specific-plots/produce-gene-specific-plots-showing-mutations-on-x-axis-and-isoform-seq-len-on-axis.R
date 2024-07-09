# Load the required library
library(ggplot2)

# Define file paths
input_file <- "enrichment_filtered_with_gene_total_variant_count.txt"
output_dir <- "gene-isoform-plots"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Extract the gene identifier (part before the hyphen) from the Identifier column
data$Gene <- sub("-.*", "", data$Identifier)

# Calculate the fraction of Variation.count / Gene_Total_Count
data$Fraction <- data$Variation.count / data$Gene_Total_Count

# Filter out isoforms with Variation.count equal to 0
data <- subset(data, Variation.count != 0)

# Get unique genes
unique_genes <- unique(data$Gene)

# Create a plot for each gene and save as PDF if it contains more than one data point
for (gene in unique_genes) {
  # Subset data for the current gene
  gene_data <- subset(data, Gene == gene)
  
  # Check if there is more than one data point
  if (nrow(gene_data) > 1) {
    # Determine which isoform is the main one (without hyphen)
    gene_data$IsMainIsoform <- ifelse(grepl("-", gene_data$Identifier), FALSE, TRUE)
    
    # Order isoforms by seq_len
    gene_data <- gene_data[order(gene_data$seq_len), ]
    
    # Create the plot
    p <- ggplot(gene_data, aes(x = reorder(Identifier, seq_len), y = Fraction)) +
      geom_point(aes(color = IsMainIsoform), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
      labs(title = paste("Dot Plot for Gene:", gene), x = "Isoform", y = "Fraction (Variation.count / Gene_Total_Count)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save the plot to a PDF file in the specified directory
    pdf(file = file.path(output_dir, paste0(gene, "_dot_plot.pdf")))
    print(p)
    dev.off()
  }
}

# Print a message indicating completion
cat("Dot plots for each gene have been saved as PDFs in the", output_dir, "directory.\n")
