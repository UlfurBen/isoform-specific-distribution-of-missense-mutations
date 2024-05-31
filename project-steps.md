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
Then I will use the total number of isoforms n and loop over n each time loading the https://rest.uniprot.org/uniprotkb/O43918-3.fasta website link with the O43918-n different each time. Then downloading this fasta file and manipulating it will be easier considering I have the available data.
My only worry is that I will only have sequence information and no pathogenicity information for further analysis.
With good documentation I will be able to review the available code and add to it to find more information.
