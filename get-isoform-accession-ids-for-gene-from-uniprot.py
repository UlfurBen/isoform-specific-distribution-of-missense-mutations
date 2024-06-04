import requests
import pandas as pd

# UniProt API endpoint for isoform information
uniprot_url = 'https://rest.uniprot.org/uniprotkb/'
kmt2d_uniprot_id = 'O14686'  # UniProt ID for KMT2D

# Function to get isoform accession numbers from UniProt
def get_isoform_accessions(uniprot_id):
    url = f"{uniprot_url}{uniprot_id}?format=json"
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Error fetching data from UniProt: {response.status_code}")
    
    data = response.json()
    isoforms = []
    
    if 'comments' in data:
        for comment in data['comments']:
            if comment['commentType'] == 'ALTERNATIVE PRODUCTS':
                for isoform in comment['isoforms']:
                    isoform_accession = isoform['isoformIds']
                    isoforms.append(isoform_accession)
    
    return isoforms

# Retrieve isoform accession numbers for KMT2D
isoform_accessions = get_isoform_accessions(kmt2d_uniprot_id)

# Print isoform accession numbers
print("Isoform accession numbers for KMT2D:")
for accession in isoform_accessions:
    print(accession)
