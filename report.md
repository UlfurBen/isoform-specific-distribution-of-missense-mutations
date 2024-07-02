# week 1
I learned ways to gather data and discovered databases that could be useful such as UniProt, ensembl and ClinVar. I tried using entrez-direct for gathering mutation data from clinvar with little progress.
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
# week 4
It is better to work without api and to download the files containing all information on uniprot.
I have created a command pipeline that can search in an isoform fasta file from uniprot for the human isoforms given a gene name.
These I can then use to search in variant .txt file (containing all variants stored on uniprot from many different databases) and can count the number of mutations of a given isoform id from each db individually or from all dbs.
I just need to finalise the pipeline and then I can try running it with the 300 genes on https://www.epigeneticmachinery.org/ website and visualize the results and analyze the post translational modification.
# week 5
I filtered the variant database to only include missense variants from ClinVar and with an rs identifier and not RCV to avoid duplicates.
I calculated the number of isoform specific variants in the database after filtering to only include ClinVar and missense variants and I that the number is around 500,000.
I counted the number of isoforms that don't have missense variant entries in the ClinVar database and the number is:
I counted and labeled the variants which have an variant associated amino acid categorical change.

# week 6
I met with Kaan to discuss the project. I added and organised files on Github and made the workflow required to get the desired results more understandable.
I see that the awk statement causes some errors in the enrichment calculations as the end file contains fewer lines (isoforms of genes) than there are genes I looked up (I also looked up isoforms but the canonical is always present, I checked).

Enrichment ratio recalculate. Count of missense variants from ClinVar for isoform of interest / (total number of missense variants in gene of interest for all isoforms of gene of interest) and multiply this number/fraction by isoform amino acid seq length. Make the graph a dot plot and order in descending order. Make dot plot for each gene showing calculation for each isoform and highlighting the canonical isoform in each dot plot for each gene. Automate this graphing process and save as pdf or svg (both can be modified in inkscape). Do calculation only for EM genes, separate dot plot as pdf or svg.
2. Count ClinVar missense mutation ranges for isoforms, adjust to show as dot plot and start with 0 and then start ranges (1-10 or 1-20 or 1-50 adjust at will). Y axis should be number of isoforms and x axis should number of missense variants (ranges start with 0) for EM genes.
3. Recalculate property change for all missense variants from ClinVar and show the same calculation as a stack bar plot displaying EM and non EM genes percentage on y axis. Show percentage on y axis (percentage is derived from number of variants of a specified categorical change of the total amount of variants).
EMs = (number of missense variants of a certain categorical change / total amount of missense variants) x 100 %.
Same for non EMs.
4. Refilter variant file to only contain missense variants from ClinVar (don’t filter the RCV away), also filter to only EM genes in appropriative calculations.
Calculate pathogenicity percentage/count of the total for EM genes.
5. On github be specific in naming for folders and files and show clearer the input and output and show workflow. Add plots to github as a pdf file. Organize github properly.
6. Check ClinVar missense variants for only EM genes (295 genes in total), number of pathogenic, benign etc. Check all missense variants in all databases (not only ClinVar but also others) for all EM genes (295 genes in total), in terms of number of pathogenic variants, benign variants etc.
