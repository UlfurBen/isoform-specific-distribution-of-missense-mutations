import os
import requests
import json
import pandas as pd

# UniProt API endpoint for isoform information
uniprot_url = 'https://rest.uniprot.org/uniprotkb/'

# download csv file from https://www.epigeneticmachinery.org/ and change csv_file_path accordingly
csv_file_path = "/Users/ulfurfjolnisson/Downloads/The Epigenetic Machinery.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Extract all rows from the "UniProt" column and save them to a list
UniProt_ids = df["UniProt"].tolist()

# Display the DataFrame
print(UniProt_ids)

# Retrieve isoform accession numbers for genes with UniProt_id and save their .# json files
for UniProt_id in UniProt_ids:
    requestURL = f"https://www.ebi.ac.uk/proteins/api/variation/{UniProt_id}?format=json"

    response = requests.get(requestURL)
    data = response.json()
    print(data)
