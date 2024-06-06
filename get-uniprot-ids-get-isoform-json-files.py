import os
import requests
import json
import pandas as pd

# UniProt API endpoint for isoform information
uniprot_url = 'https://rest.uniprot.org/uniprotkb/'

# Path to your CSV file
csv_file_path = "/Users/ulfurfjolnisson/Downloads/The Epigenetic Machinery.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Extract all rows from the "UniProt" column and save them to a list
UniProt_ids = df["UniProt"].tolist()

# Display the DataFrame
print(UniProt_ids)

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
for UniProt_id in UniProt_ids:
    isoform_accessions = get_isoform_accessions(UniProt_id)
    # Strip the quotes from isoform accessions
    stripped_isoform_accessions = [accession[0].strip("'") for accession in isoform_accessions]

    print(stripped_isoform_accessions)
    # Iterate over stripped isoform accessions and fetch data
    for stripped_isoform_accession in stripped_isoform_accessions:
        url = f"https://www.ebi.ac.uk/proteins/api/variation/{stripped_isoform_accession}?format=json"
        response = requests.get(url)
        data = response.json()
        print(data)
        
        folder_name = UniProt_id
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        # Save the data to a file named after the stripped_isoform_accession variable value
        file_name = os.path.join(folder_name, f"{stripped_isoform_accession}.json")
        with open(file_name, 'w') as f:
            json.dump(data, f)

# # Create a new directory named after gene name if it doesn't exist
# folder_name = 'kmt2d'
# if not os.path.exists(folder_name):
#     os.makedirs(folder_name)

