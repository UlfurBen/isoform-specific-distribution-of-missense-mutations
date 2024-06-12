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

# 06.06.24

I downloaded .csv file from https://www.epigeneticmachinery.org/ site and loaded UniProt IDs as a list. I iterated over the UniProt IDs and found isoform accession IDs for their isoforms. I got and downloaded isoform .json files using the isoform accession IDs. Code for this is in this repository on Github. I now have the code to fetch variant and isoform information from UniProt.
Next step is to learn GREP to get the part of the data I am interested in using.
I also have to know which information will be of use for the enrichment analysis.
I must add that this part felt hardest considering I was using novel and not documented ways at first but ended up using well documented ways with the api.

# 07.06.24

I now have the code required to download variation information from list of gene UniProt IDs and save the result to files.
I also have corresponding file structures that downloads isoform information from the accession IDs of the isoforms which are fetched in the same script from the UniProt gene names.
I tried using grep on .json file but realized I need jq tool to separate json file data from one line to many lines so grep can accurately analyze it.
I wanted to incorporate ensembl information about the isoforms of the genes into my dataset.
I now have a script in this repo that gets ensembl information with the api from ensembl IDs which I got from the gene names from the epigenetic-machinery.org site. I had to write a script to find all the ensemble ids from the gene names which I then added to the .csv file I downloaded from the epigenetic-maginery.org site.

This weekend and next week will focus on getting mutation info from clinvar which has proven difficult because of their hard to use api (I have to use biopython to use entrez to get the information and the bugs are many as of now).
I will mostly focus next week on fetching the information in the files to make clean files for enrichment analysis.

I am pretty ahead of schedule and maybe I will use novel ways later to answer the research question: "are there isoforms/isoform-areas with unusually high pathogenic mutations.

# 10.06.24

Steps:
1. Get data , uniprot contains mutation data, as well as isoform data but not the start and end data of the isoforms only the length of the isoforms. I have ensembl data which contains exon info for a given gene with chromosome, length, start and end data. I have working code to retrieve this data given a list of gene names and database IDs.
2. Next I'll need to filter the relevant components of the data for plotting. The idea is to show on the Y-axis the number of pathogenic or highly-pathogenic missense mutations and on the x-axis show the isoform and/or the exon where the enrichment is located.

# 11.06.24

Now I have the data I want to filter: isoforms ids and sequence data from isoform files. Exon data from ensembl files.
I also want to filter mutation type, genomic location and frequency (population frequency) from variation data.
In both of these filters I want to save the resulting data to .csv files.
I should propably have gene names and/or isoform accession ids present in all files so that I can align the data correctly.
It is proving to be a headache to fetch the data from clinvar since the ftp site is not useful/accessible/thorough enough and the entrez direct cannot fetch the correct data for a given gene, not to mention the waiting time needed when quering for data which causes you to make many batches or requests instead of one long one.
I therefore decided to use myvariant.info instead and at first glance th data is thorough and accurate and easy to get.

# 12.06.24

Kaan explains lof intolerance (pli)
Manually check exon info (get the info manually to compare to automatic retrieval results)

*Understand json file structure*

*Isoform specific variant info*

Understand pli score
Vus (variant of uncertain signifance)
