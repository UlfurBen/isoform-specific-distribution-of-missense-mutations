# week 1
I learned ways to gather data and databases that would be useful. I tried using entrez-direct for gathering mutation data from clinvar with little progress.
I then tried using ucsc genome browser online where I uploaded some tracks of mutations and visualised them but realised that the work need to download data individually
and then to uploaded to the browser would prove inefficient and cumbersome. I therefore learned to use APIs instead and started using the REST API which is well compatable
with UniProt databases.
# week 2
I tried using chatgpt to write python code that followed steps I defined towards finding the right data. I tried doing the same with clinvar, using different python packages
but the data wasn't better than on UniProt. UniProt contains some clinvar data which is extremely useful.
I started making code that gets the ensembl ids using the gene names on epigenetic machinery .org website and then to get exon data using the ids from ensembl using the API.
# week 3
I finalised and stitched together bits of scripts to use the csv file with a list of 300 genes we want to analyse and lookup and download relevant json files from UniProt using the API.
I downloaded the variant information of the genes, the isoforms of the genes and the exons of the genes and saved the information within the json files to csv files for further analysis and plotting.
I managed to get some code working in python that takes exon and variant csv files for some genes (only used 4 in this analysis for quickness but it's easy to scale) and plotted the enrichment of missense
mutations across the exons.
Kaan wrote a question I sent to UniProt staff about getting isoform specific variants for isoforms of genes that aren't computational isoforms since this is the wall I'm facing in finalising this project.
