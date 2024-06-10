import pandas as pd
import requests
import os
import json

# Path to your CSV file
csv_file_path = "/Users/ulfurfjolnisson/Downloads/The Epigenetic Machinery.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Extract all rows from the "Ensembl ID" column and save them to a list
ensembl_ids = df["Ensembl ID"].tolist()

# Extract all rows from the "Gene Name" column and save them to a list
gene_names = df["Gene Name"].tolist()

# Ensure the output directory exists
output_dir = "isoform_information"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Ensembl REST API base URL
ensembl_base_url = "https://rest.ensembl.org"

# Function to get isoform information
def get_isoform_info(ensembl_id):
    url = f"{ensembl_base_url}/lookup/id/{ensembl_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error fetching data for {ensembl_id}: {response.status_code}")
        return None

# Iterate over each Ensembl ID and Gene Name
for ensembl_id, gene_name in zip(ensembl_ids, gene_names):
    # Get isoform information
    isoform_info = get_isoform_info(ensembl_id)
    
    if isoform_info:
        # Create output file path
        output_file_path = os.path.join(output_dir, f"{gene_name}_ensembl.json")
        
        # Write isoform information to the file
        with open(output_file_path, 'w') as file:
            json.dump(isoform_info, file, indent=4)
        
        print(f"Saved isoform information for {gene_name} to {output_file_path}")

print("Isoform information fetching complete.")
