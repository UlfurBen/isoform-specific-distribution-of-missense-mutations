# 28.05.24

Proposed methods for finding missense mutation enriched areas across epigenetic isoforms:
Integrative genome viewer
Flibase

Ensembl

Uniprot

Ucsc

epigenetic machinery

Snpeff

Rstudio will be used for data gathering, analysis and visualization.

Tutorials used:
	https://brouwern.github.io/lbrb/introduction-to-biological-sequences-databases.html#the-fasta-file-format


Fasta data cleaning:

	fasta_cleaner() under the package compbio4all from https://brouwern.github.io/lbrb/worked-example-building-a-phylogeny-in-r.html#download-necessary-packages 


# 29.05.24

Use ClinVar api to fetch and download gene variant data of epigenetic gene hosted by epigeneticmachinery.org

Use the ftp clinvar website to get data from the xml format. This is the website: ​​https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ 

Use web scraping to find the variant information from clinvar.

Ensembl vep

Command to get the variant ids:
esearch -db clinvar -query "missense_variant AND KMT2D["Gene"]" | efetch -format xml > GENE_NAME_missense_variants.xml

# 30.05.24

Fetch out all isoform sequences for kmt2d + computational ones as json format

Check for pairwise alignment between api and manual lookup on uniprot

What is downstream process of json file format

Use perl: perl executable.pl File.json > File.txt

Mark down file format (learn) hosted on github with all steps and intermediary results documented

used this link:
https://www.uniprot.org/help/api_queries 
in order to get an idea about structuring uniprot api_queries to get isoform sequences

# 31.05.24

The AI generated code doesn't seem to work as intended, never adjusting itself to produce code that finds the isoform sequence information.
My idea is to use python web scraping using the BeautifulSoup library to find the element on the uniprot entry site that contains information about the total amound of isoforms.
Then I will use the total number of isoforms n and loop over n each time loading the https://rest.uniprot.org/uniprotkb/O43918-3.fasta website link with the O43918-i for i = 0; i < n; i++. Then downloading this fasta file and manipulating it will be easier considering I have the available data.
My only worry is that I will only have sequence information and no pathogenicity information for further analysis.
With good documentation I will be able to review the available code and edit it to find further information.
Automating this process is difficult because of my unfamiliarity with python web scraping. Therefore I will try to manually download the missense mutation and isoform data for the aire gene and analyze that using R.
Hopefully with this knowledge and understanding I will be able to make a program in python that accesses this data simultaneously for all known genes.

# 01.06.24

I learned how to edit a data frame in rstudio. The problem I got when reading the data from the mutation data frame arose because of N/A error.

# 03.06.24

I sent a message to uniprotkb staff with a question regarding the use of their api.

ISOformSwitchAnalyzeR package could be useful for visualizing the distribution of missense mutations across isoforms.

Hans proposed the idea to only plot the enrichment of mutations and confine ourselves with the areas with the highest enrichment.

Another idea is to programmatically erase all isoform areas from the data that overlaps with all other isoforms and confine ourselves with the areas that are present only in few isoforms.

With the enrichment information we can limit our results to only areas with a mutation density above a certain value.

I wrote made ChatGPT write R code to save first column information from ClinVar mutations data frame.

# 04.06.24

I add .py files to github. One file contains code to get accession ids for gene of interest. Another file uses accession ids to get start, end and length information. Other files are in R also. They do simple tasks like manipulate mutation data from ClinVar by saving first and third column to new data file.

# 05.06.24

Reverse engineering this project is easier to understand. Final step consists of plotting in R using ggplot2 the distribution of mutation on isoforms. Data is structured with column 1 containing isoform ID (x axis) and column 2 containing mutation enrichment (y axis).
