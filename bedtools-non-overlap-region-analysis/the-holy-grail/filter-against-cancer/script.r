# List of gene names to filter
gene_list <- c("AEBP1", "C1orf189", "APTX", "SLC25A22", "QTRT1", "SNX15", "LZIC", 
               "EIF5B", "NFIB", "LEPREL2", "SLC9A3R2", "SCGB1D1", "SLC12A1", 
               "ABCG4", "LRRC37A", "ALX3", "EIF4G3", "GALT", "ITGAE", "NFE2L3", 
               "BICD2", "RAB4A", "LPCAT3", "RABEP2", "NEURL4", "KCNH3", "PKMYT1", 
               "PLEKHH3", "CACNG2", "SH3BP2", "XPO5", "NPRL3", "TPM2", "DOHH", 
               "MST1", "B9D2", "LSM11", "PFKFB4", "TOR2A", "IL27RA", "SLC6A5", 
               "GRHPR", "TTC16", "ZER1", "SPC24", "MMP8", "SMCO1", "UBR1", 
               "CCNA2", "ZSWIM5", "CYP2A6", "DKFZP434E1119", "LLGL1", "C6orf222", 
               "MTMR14", "ANKH", "TRIOBP", "KLK2", "DDIT4L", "PBX3", "CAMKMT", 
               "LRP4", "MC2R", "MRPS10", "FANCG", "CENPQ", "GIPC3", "PIEZO1", 
               "RP11-894J14.5", "TEK", "GNA11", "NPR2", "LRP4", "KCNV2", "LSM7", 
               "TDRD12", "MLLT1", "PSKH1", "ITGAE", "TUBB2B", "SURF6", "TUBB2A", 
               "FAM105B", "FLII", "FOXC1", "IFITM5", "LRRC56", "PEX7", "HCN1", 
               "EEF2", "PIEZO1", "NHLRC1", "LAMA1", "UTP20", "RUSC1", "HOOK2", 
               "DPCR1", "SLC1A1", "OXCT1", "AC138655.1", "TYRP1", "POLRMT", 
               "SIPA1L3", "FGF10", "WIZ", "PALB2", "RP11-433C9.2", "KISS1R", 
               "KRT81", "ZNF408", "AC174470.1", "CSH1", "GCM2", "TMEM151A", 
               "DMRT1", "PLAA", "FBXW9", "PTGER2", "AC132872.2", "PINX1", 
               "AMH", "HSPA1A", "OR1F1", "ADAT3", "MAFA", "BMP6", "C1orf233", 
               "FFAR1", "DMRT3", "PLEKHO2", "GPR179", "AP5B1", "KIAA1161", 
               "MUC12", "FAM230A")

# Read the BED file
formatted_variants <- read.table("formatted_variants_likely_p_p.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Add column names
colnames(formatted_variants) <- c("chr", "chr_start", "chr_end", "transcript_id", "gene_name", "mut_count", "enrichment_ratio")

# Filter the dataframe based on the gene list
filtered_variants <- formatted_variants[formatted_variants$gene_name %in% gene_list, ]

# Write the filtered dataframe back to the BED file
write.table(filtered_variants, "formatted_variants_likely_p_p.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
