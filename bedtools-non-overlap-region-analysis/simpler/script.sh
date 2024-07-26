# Step 1: Generate genome file
awk '{if($3 > max[$1]) {max[$1] = $3}} END {for (chr in max) print chr, max[chr]}' Homo_sapiens.GRCh37.87_without_headers.bed > genome_file.txt

# Step 2: Find non-overlapping regions
bedtools complement -i Homo_sapiens.GRCh37.87_without_headers.bed -g genome_file.txt > non_overlapping_regions.bed

# Step 3: Intersect with variants and count
bedtools intersect -a non_overlapping_regions.bed -b homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed -c > intersected_counts.bed

# Step 4: Calculate enrichment
awk '{
    length = $3 - $2;
    enrichment = $4 / length;
    print $0 "\t" length "\t" enrichment
}' intersected_counts.bed > intersected_counts_with_enrichment.bed
