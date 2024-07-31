# Load the required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# Define file paths
# Input file is the name on my local computer but on elja the file's name is enrichment_filtered_with_gene_total_variant_count.txt
input_file <- "enrichment.txt"
epigenetic_file <- "The-Epigenetic-Machinery.csv"
output_dir <- "gene-isoform-plots"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")
epigenetic_data <- read.csv(epigenetic_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Calculate the fraction of Variation.count / Gene_Total_Count
data$Fraction <- data$Variation.count / data$Gene_Total_Count * 100

# Get unique gene names
unique_gene_names <- unique(data$Gene_Name)

# Create a plot for each gene and save as PDF
for (gene_name in unique_gene_names) {
  # Subset data for the current gene name
  gene_data <- subset(data, Gene_Name == gene_name)
  
  # Create a directory for the current gene name if it does not exist
  gene_dir <- file.path(output_dir, gene_name)
  if (!dir.exists(gene_dir)) {
    dir.create(gene_dir)
  }
    
  # Filter out data points where either seq_len or Fraction is 0
  # gene_data <- gene_data[gene_data$seq_len != 0 & gene_data$Fraction != 0, ]
    
  # Determine which isoform is the main one (without hyphen)
  canonical_identifiers <- epigenetic_data$UniProt
  gene_data$IsMainIsoform <- gene_data$Identifier %in% canonical_identifiers

  
  # Order isoforms by seq_len
  gene_data <- gene_data[order(gene_data$seq_len), ]
    
  # Create the dot plot
  p_dot <- ggplot(gene_data, aes(x = seq_len, y = Fraction)) +
    geom_point(aes(color = IsMainIsoform), size = 3) +
    geom_text_repel(aes(label = Identifier), size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
    labs(title = paste("Dot Plot for Gene:", gene_name), x = "Isoform Sequence Length", y = "Fraction (Variation.count / Gene_Total_Count) * 100") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  # Save the dot plot to a PDF file in the specified directory for the gene name
  pdf(file = file.path(gene_dir, paste0(gene_name, "_dot_plot.pdf")))
  print(p_dot)
  dev.off()
    
  # Prepare data for the new stacked bar plot
  gene_data_new_long <- gene_data %>%
    mutate(Total_Count = Benign + Likely.benign + VUS + Likely.pathogenic + Pathogenic) %>%
    pivot_longer(cols = c(Benign, Likely.benign, VUS, Likely.pathogenic, Pathogenic), 
                 names_to = "VariantType", values_to = "Count") %>%
    group_by(Gene_Name, Identifier, VariantType, Total_Count) %>%
    summarize(Count = sum(Count)) %>%
    ungroup() %>%
    mutate(Percentage = (Count / Total_Count) * 100)
    
  # Check if there is data for the new stacked bar plot
  if (nrow(gene_data_new_long) > 0) {
    # Create the new stacked bar plot
    p_new_bar <- ggplot(gene_data_new_long, aes(x = Identifier, y = Percentage, fill = VariantType)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Variant Type Distribution for Gene:", gene_name), x = "Identifier", y = "Percentage of Total Variants") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    # Save the new stacked bar plot to a PDF file in the specified directory for the gene name
    pdf(file = file.path(gene_dir, paste0(gene_name, "_variant_type_distribution.pdf")))
    print(p_new_bar)
    dev.off()
  }
}

# Print a message indicating completion
cat("Dot plots and variant type distribution plots for each gene have been saved as PDFs in the", output_dir, "directory.\n")
