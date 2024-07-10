# Load the required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# Define file paths
input_file <- "enrichment.txt"
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
    # Create a directory for the current gene if it does not exist
    gene_dir <- file.path(output_dir, gene)
    if (!dir.exists(gene_dir)) {
      dir.create(gene_dir)
    }
    
    # Filter out data points where either seq_len or Fraction is 0
    gene_data <- gene_data[gene_data$seq_len != 0 & gene_data$Fraction != 0, ]
    
    # Determine which isoform is the main one (without hyphen)
    gene_data$IsMainIsoform <- ifelse(grepl("-", gene_data$Identifier), FALSE, TRUE)
    
    # Order isoforms by seq_len
    gene_data <- gene_data[order(gene_data$seq_len), ]
    
    # Check again if there is more than one data point after filtering
    if (nrow(gene_data) > 1) {
      # Create the dot plot
      p_dot <- ggplot(gene_data, aes(x = seq_len, y = Fraction)) +
        geom_point(aes(color = IsMainIsoform), size = 3) +
        geom_text_repel(aes(label = Identifier), size = 3) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
        labs(title = paste("Dot Plot for Gene:", gene), x = "Isoform Sequence Length", y = "Fraction (Variation.count / Gene_Total_Count)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Save the dot plot to a PDF file in the specified directory for the gene
      pdf(file = file.path(gene_dir, paste0(gene, "_dot_plot.pdf")))
      print(p_dot)
      dev.off()
    }
    
    # Prepare data for stacked bar plot, filtering out isoforms with any 0 values in specified columns
    gene_data_long <- gene_data %>%
      filter(Benign != 0, Likely.benign != 0, VUS != 0, Likely.pathogenic != 0, Pathogenic != 0) %>%
      select(Gene, Identifier, Gene_Total_Count, Benign, Likely.benign, VUS, Likely.pathogenic, Pathogenic) %>%
      pivot_longer(cols = c(Benign, Likely.benign, VUS, Likely.pathogenic, Pathogenic), names_to = "VariantType", values_to = "Count") %>%
      group_by(Gene, Identifier, Gene_Total_Count, VariantType) %>%
      summarize(Count = sum(Count)) %>%
      ungroup() %>%
      mutate(Percentage = (Count / Gene_Total_Count) * 100)
    
    # Check if there is data for the stacked bar plot
    if (nrow(gene_data_long) > 0) {
      # Create the stacked bar plot
      p_bar <- ggplot(gene_data_long, aes(x = Identifier, y = Percentage, fill = VariantType)) +
        geom_bar(stat = "identity") +
        labs(title = paste("Variant Distribution for Gene:", gene), x = "Identifier", y = "Percentage of Gene_Total_Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Save the stacked bar plot to a PDF file in the specified directory for the gene
      pdf(file = file.path(gene_dir, paste0(gene, "_stacked_bar_plot.pdf")))
      print(p_bar)
      dev.off()
    }
  }
}

# Print a message indicating completion
cat("Dot plots and stacked bar plots for each gene have been saved as PDFs in the", output_dir, "directory.\n")
