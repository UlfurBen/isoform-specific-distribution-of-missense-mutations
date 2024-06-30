#!/bin/bash

# Define the input and output file paths
file1="cleaned_all_300_calculate_missense_variant_enrichment_within_isoforms_with_lengths_for_first_isoform_of_genes.txt"
file2="cleaned_all_300_calculate_missense_variant_enrichment_within_isoforms_with_lengths_for_rest_of_isoforms_of_genes.txt"
output_file="merged_calculate_missense_variant_enrichment_within_isoforms_with_lengths_all_300.txt"

# Merge the two files into the output file
cat "$file1" "$file2" > "$output_file"

# Print a message indicating completion
echo "Files have been merged and saved to $output_file."

