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

# 14.06

Use ENST, ENSG abd ENSP for cross referencing variations and isoforms.
Make bar plot using ftp file to show the enrichment of mutations across isoforms but also show the db source for percentage of mutations.
Separate plot using only clinvar and showing percentage of benign to pathogenic (lb, b, lp, p and vus).
Use rgraphgallery for bar plot.
Bar plot showing on x axis the db source and on y axis the total variant number.

We use
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/
for downloading variants and use unix commands on the terminal to look for data we want.

I can use isoform fasta file of all isoforms and search for isoforms of a gene by gene name. Then I can find isoform id which is present in variant header in variant file. I can search for variant with header including isoform id and count number of lines since each line contains new variant. Number of lines returned for isoform id equals number of variants per that isoform.
I can then do this for all gene names, search for their name in isoform fasta file and find their isoform id. Then search in variant txt file for variant entries with their isoform id in header.

Only problem I realise with this approach is if isoform fasta file includes amino acid sequence with subsequence matching gene name, for example AIRE which after grep-ing I found to be included in several amino acid sequences of isoforms.

# 18. june

I solved my proposed "problem", in the fasta zip file I always search for {gene_name}_HUMAN so the suffix would always be _HUMAN.

I visited 
https://www.uniprot.org/help/downloads
and downloaded Isoforms sequences as fasta file.
I then ran 

gunzip -c uniprot_sprot.fasta.gz | grep KMT2A_HUMAN

which returned
>sp|Q03164|KMT2A_HUMAN Histone-lysine N-methyltransferase 2A OS=Homo sapiens OX=9606 GN=KMT2A PE=1 SV=5
of which Q03164 is the input for the next command.

From the link I provided above I ran

gunzip -c homo_sapiens_variation.txt.gz | grep O43918 | wc -l

which returns the mutations associated with the O43918 gene and the number of mutations (which equals then number of lines of output.

The file "humsavar.txt" contains the variations for the first isoform of each gene.
The file "uniprot_sprot_varsplic.fasta.gz" contains isoform ids of each gene except for the first one.
The file "homo_sapiens_variation.txt.gz" contains the variations for the rest of the isoforms (all but the first one).

Workflow example for KMT2A:
touch output.txt
echo -e "\nQ03164-2: variation number:" >> output.txt
gunzip -c uniprot_sprot_varsplic.fasta.gz | grep KMT2A_HUMAN
gives me "Q03164-2" and "Q03164-3"
echo -e "\nQ03164-2: variation number:" >> output.txt
gunzip -c homo_sapiens_variation.txt.gz | grep Q03164-2 | wc -l >> output.txt
echo -e "\nQ03164-3: variation number:" >> output.txt
gunzip -c homo_sapiens_variation.txt.gz | grep Q03164-3 | wc -l >> output.txt
cat output.txt
To get the canonical isoform I run:
gunzip -c homo_sapiens_variation.txt.gz | grep 'Q03164' | grep -v 'Q03164-' | wc -l
For KMT2A_HUMAN isoforms the output file I would graph with is:
Name,Count
Q03164-1,6259
Q03164-2,5939
Q03164-3,6864
To get the gene id for running in variation file I run:

gunzip -c uniprot_sprot_varsplic.fasta.gz | grep KMT2A_HUMAN


and I get

>sp|Q03164-2|KMT2A_HUMAN Isoform 2 of Histone-lysine N-methyltransferase 2A OS=Homo sapiens OX=9606 GN=KMT2A
>sp|Q03164-3|KMT2A_HUMAN Isoform 3 of Histone-lysine N-methyltransferase 2A OS=Homo sapiens OX=9606 GN=KMT2A



_Full guide_:


gunzip -c uniprot_sprot_varsplic.fasta.gz | grep 'KMT2A_HUMAN' > temp.txt

awk -F '[|-]' '{print $2}' temp.txt | sort -u > temp_identifiers.txt

while read -r identifier; do
  count=$(gunzip -c homo_sapiens_variation.txt.gz | grep -w "$identifier" | wc -l)
  echo "${identifier},${count}" >> output-3.txt
done < temp_identifiers.txt

awk -F '|' '{print $2}' temp.txt | sort -u > temp_identifiers_2.txt

while read -r identifier; do
  count=$(gunzip -c homo_sapiens_variation.txt.gz | grep -w "$identifier" | wc -l)
  echo "${identifier},${count}" >> output-3.txt
done < temp_identifiers_2.txt

rm temp.txt temp_identifiers.txt temp_identifiers_2.txt

# 19. june

I can create the pipeline entirely in R packaging the shell commands inside the system() function.
I intend to use gene names and append _HUMAN to them and then to use those edited names to grep the isoform fasta file for isoform ids and to use those as identifiers to grep the variant .txt file and count the mutation number per isoform.
I used elja for the first time and submitted an sbatch job that calculates the number of mutations for the isoforms of 6 genes (just used 2 cpu node, 4 gb of ram and 1 hour max run time.

# 20. june

The results from the sbatch job were as expected.
I decided to create another R file. That file filters the data from the first R file so that it subtracts from isoform one count the count of the other isoforms of the same gene. This was necessary because I double counted in the first R script.
I decided to edit the first R script (iterative-loop.R) to calculate accurately the counts of the first isoform.
If this works then I can scale it to also calculate the counts for pathogenicity and percentage from each of the 10 database sources.
I managed to fix the script and now I can add more logic to add more information to the output file that I can use to plot graphs.

# 21. june
For next week I will present the source file contents for Kaan at 10 am Monday.
I want to:
	count the total amount of variants in homo_sapiens_variation.txt.gz
 	count the number of unique isoforms in homo_sapiens_variation.txt.gz
 	count the number of variants per database (10 databases in total)
  	calculate the percentage of variants per database per isoform from the total isoform variant number
   	count the number of amino acid changes that have polarity (after mutation) and also have high pathogenicity, 		do this per isoform specific variant


Results:
	There are 52 985 969 total variants in homo_sapiens_variation.txt.gz
 
 	There are 95 368 unique splice isoforms in homo_sapiens_variation.txt.gz
  
  	The number of missense variants per database:
   	ClinVar,   dbSNP,     ESP,        ExAC,      TOPMed,     gnomAD,    NCI-TCGA Cosmic,cosmic curated, 1000Genomes
   	8 395 578, 4 181 656, 3 914 220, 17 890 491, 23 432 606, 32 588 009, 1 576 481,      1 225 408,     3 653 273

# 24. june
In the meeting with Kaan we talked about looking at GTEX database for tissue specific expression information. Is this splice isoform only expressed in the brain? Is it only expressed in the liver? For instance.

✅ We talked about saving to a new file (called homo_sapiens_variation_missense_ClinVar.txt) only the entries containing only missense variants that are from the ClinVar database from the original variant file from downloaded from uniprot.

✅ We talked about taking homo_sapiens_variation_missense_ClinVar.txt and counting the number of isoform specific variants using the source DB ID (which starts with rs).

We talked about counting in homo_sapiens_variation_missense_ClinVar.txt the number of isoforms with no variants associated with them but to find all documented isoforms we use uniprot_sprot_varsplic.fasta.gz.

We talked about how showing numbers in a plot is more insightful and expressive. I will do that from now on as well as create a google slides presentation for clarity.


- I created numerous files that contain filtered the variant data according to database, variant type, RCV vs rs id etc.
- I created on elja the file "trial_count_polarity_or_charged_change.R" that correctly checks whether a variant results in a polar amino acid.

# 25. june
I have some results in the following files stored on elja:
	count_of_isoforms_without_variant_entry.txt
 	trial_count_polarity_or_charged_change.txt

# 26. june
I'm still working on the missense variant amino acid property change R script that checks the property change of the original a.a. vs the new a.a. in a missense variant.
I got the property change script to work.

✅ trial_count_polarity_or_charged_change_property_change_only.R
✅ Amino_acid_change_with_property_change_updated.txt
   use uniprot_sprot_varsplic.fasta.gz to search calculate_missense_variant_enrichment_within_isoforms.txt for 
   isoform ids and calculate and save the isoform sequence length

   In calculate_missense_variant_enrichment_within_isoforms_with_lengths.txt I have isoform 1 of each gene sequence length.
   In calculate_missense_variant_enrichment_within_isoforms_with_lengths.txt I have the sequence length of the rest of the isoforms.
✅ In "merged_calculate_missense_variant_enrichment_within_isoforms_with_lengths.txt" I have all isoforms with their length and ClinVar missense variation count.
✅ In ~/Downloads/mutation_sequence_length_fraction.R on my computer I have plot to show the mutation_count to sequence_length ratio for each isoform
I have to change the ratio to mutation_count/(sequence_length x gene_isoform_number)

For tomorrow:
1. change ratio to include the total gene isoform number
2. Count the average isoform number for EM genes (all and only EM genes are in "homo_sapiens_variation_missense_ClinVar_Reference_SNP_EM_genes.txt"
3. Do all calculations across all 300 EM genes (I only did the analysis across 6 EM genes)
4. Count number of phenotype information (pathogenic, vus etc.) in EM variant file.
5. Count number of isoforms which have variant count in a certain range (like 5 to 10, 10 to 20 etc.).
