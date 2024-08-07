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
I have solved the enrichment analysis problem where the program would miss calculating some isoforms. Now I have a stacked bar plot that show benign, vus, likely benign, pathogenic and likely pathogenic in each gene and shows the bars in descending order dependent on the total variant counts.

# week 7
I had an idea to find the area length of the aggregation of the missense mutations within each but the length was much bigger than the isoform length itself...
I learned to work with files in R instead of using grep, using fread tends to be quicker.
I learn to organize code better on github and I make gene dependent isoform dot plots showing on y axis the isoform variant counts divided by the total gene variant counts and on x axis the isoform sequence length. I will also make stacked bar plot showing p, lp, vus, lb and b pathogenicity information.
