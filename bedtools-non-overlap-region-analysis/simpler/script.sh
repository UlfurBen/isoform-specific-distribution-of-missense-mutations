# Step 1: Sort and filter the BED file for standard chromosomes
awk '$1 ~ /^([0-9]+|X|Y|MT)$/' Homo_sapiens.GRCh37.87_without_headers.bed | \
    sort -k1,1V -k2,2n > Homo_sapiens.GRCh37.87_without_headers_filtered_sorted.bed

# Step 2: Generate genome file from the sorted and filtered BED file
awk '{if($3 > max[$1]) {max[$1] = $3}} END {for (chr in max) printf "%s\t%.0f\n", chr, max[chr]}' \
    Homo_sapiens.GRCh37.87_without_headers_filtered_sorted.bed > genome_file.txt

# Ensure genome file is correctly sorted
sort -k1,1V genome_file.txt > genome_file_sorted.txt
mv genome_file_sorted.txt genome_file.txt

# Verify genome file format
cat genome_file.txt

# Step 3: Find non-overlapping regions
bedtools complement -i Homo_sapiens.GRCh37.87_without_headers_filtered_sorted.bed -g genome_file.txt > non_overlapping_regions.bed

# Step 4: Reformat and sort the variant file to avoid scientific notation
awk '{printf "%s\t%.0f\t%.0f\n", $1, $2, $3}' homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed | \
    sort -k1,1V -k2,2n > homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_fixed_sorted.bed

# Step 5: Identify and remove invalid records from the variant file
awk 'NF == 3' homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_fixed_sorted.bed > \
    homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_fixed_sorted_cleaned.bed

# Step 6: Intersect with variants and count
bedtools intersect -a non_overlapping_regions.bed -b homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers_fixed_sorted_cleaned.bed -c > intersected_counts.bed

# Create the awk script file
cat << 'EOF' > calculate_enrichment.awk
{
    length = $3 - $2;
    enrichment = $4 / length;
    print $0 "\t" length "\t" enrichment;
}
EOF

# Step 7: Calculate enrichment using the awk script
awk -f calculate_enrichment.awk intersected_counts.bed > intersected_counts_with_enrichment.bed
