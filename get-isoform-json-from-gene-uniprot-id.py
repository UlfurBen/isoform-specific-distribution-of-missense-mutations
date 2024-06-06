import os
import requests
import json

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

# Strip the quotes from isoform accessions
stripped_isoform_accessions = [accession[0].strip("'") for accession in isoform_accessions]

# Create a new directory named 'kmt2d' if it doesn't exist
folder_name = 'kmt2d'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# Iterate over stripped isoform accessions and fetch data
for stripped_isoform_accession in stripped_isoform_accessions:
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{stripped_isoform_accession}?format=json"
    response = requests.get(url)
    data = response.json()
    
    # Save the data to a file named after the stripped_isoform_accession variable value
    file_name = os.path.join(folder_name, f"{stripped_isoform_accession}.json")
    with open(file_name, 'w') as f:
        json.dump(data, f)
