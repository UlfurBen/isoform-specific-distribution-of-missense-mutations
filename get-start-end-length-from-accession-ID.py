import requests
import pandas as pd

# List of isoform accession numbers for KMT2D
accession_numbers_kmt2d = ['O14686-1', 'O14686-3']

# UniProt REST API URL
uniprot_url = 'https://rest.uniprot.org/uniprotkb/'

# Function to get isoform details from UniProt
def get_isoform_details(accession_numbers):
    isoforms_data = []
    
    for accession in accession_numbers:
        url = f"{uniprot_url}{accession}?format=json"
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Error fetching data for {accession}: {response.status_code}")
            continue
        
        data = response.json()
        sequence = data.get('sequence', {}).get('value', '')
        isoform_data = {
            'accession': accession,
            'start': 1,
            'end': len(sequence),
            'length': len(sequence)
        }
        isoforms_data.append(isoform_data)
    
    return pd.DataFrame(isoforms_data)

# Retrieve isoform details
isoform_details_df = get_isoform_details(accession_numbers_kmt2d)

# Save the isoform details to a CSV file for further use
# isoform_details_df.to_csv('isoform_details_kmt2d.csv', index=False)

# Print isoform details
print(isoform_details_df)
