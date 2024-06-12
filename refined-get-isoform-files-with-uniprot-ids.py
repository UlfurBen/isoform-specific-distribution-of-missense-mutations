import os
import requests
import json
import pandas as pd

# Path to your CSV file (downloaded from the epigenetic machinery .org website
csv_file_path = "/Users/ulfurfjolnisson/Downloads/The Epigenetic Machinery.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Extract all rows from the "UniProt" column and save them to a list
UniProt_ids = df["UniProt"].tolist()

# Display the DataFrame
print(UniProt_ids)

for UniProt_id in UniProt_ids:
  # Iterate over stripped isoform accessions and fetch data
  requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{UniProt_id}/isoforms.json"

  response = requests.get(requestURL)
  data = response.json()
  print(data)

  file_name = f"{UniProt_id}_isoform_refined.json"
  # Save the data to a file named after UniProt ID
  with open(file_name, 'w') as file:
        json.dump(data, file, indent = 4)
