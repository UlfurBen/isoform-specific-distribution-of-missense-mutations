bedtools complement -i Homo_sapiens.GRCh37.87_without_headers.bed -g genome_file.txt > non_overlapping_regions.bed
bedtools intersect -a non_overlapping_regions.bed -b homo_sapiens_variation_missense_ClinVar_filtered_relevancy_no_headers.bed -c > intersected_counts.bed
awk '{
    length = $3 - $2;
    enrichment = $4 / length;
    print $0 "\t" length "\t" enrichment
}' intersected_counts.bed > intersected_counts_with_enrichment.bed
