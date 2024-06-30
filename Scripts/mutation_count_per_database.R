# Define the output file
output_file <- "mutation-count-per-database.txt"

# Count occurrences for each database
count_ClinVar <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'ClinVar' | wc -l"), intern = TRUE))
count_dbSNP <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'dbSNP' | wc -l"), intern = TRUE))
count_ESP <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'ESP' | wc -l"), intern = TRUE))
count_ExAC <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'ExAC' | wc -l"), intern = TRUE))
count_TOPMed <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'TOPMed' | wc -l"), intern = TRUE))
count_gnomAD <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'gnomAD' | wc -l"), intern = TRUE))
count_NCI_TCGA_Cosmic <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'NCI-TCGA Cosmic' | wc -l"), intern = TRUE))
count_cosmic_curated <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep 'cosmic curated' | wc -l"), intern = TRUE))
count_1000Genomes <- as.numeric(system(paste0("gunzip -c homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt | grep '1000Genomes' | wc -l"), intern = TRUE))

# Write the counts to the output file
write(paste0(count_ClinVar, ",", count_dbSNP, ",", count_ESP, ",", count_ExAC, ",", count_TOPMed, ",", count_gnomAD, ",", count_NCI_TCGA_Cosmic, ",", count_cosmic_curated, ",", count_1000Genomes),
      file = output_file,
      append = TRUE)