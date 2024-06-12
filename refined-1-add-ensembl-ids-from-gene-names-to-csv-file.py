import os
import requests
import json
import pandas as pd

csv_file_path = "/Users/ulfurfjolnisson/Documents/The Epigenetic Machinery.csv"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Extract all rows from the "Gene Name" column and save them to a list
gene_names = df["Gene Name"].tolist()

def get_ensembl_id(gene_names):
    server = "https://rest.ensembl.org"
    endpoint = "/lookup/symbol/homo_sapiens/"
    headers = {"Content-Type": "application/json"}

    gene_ids = {}
    for gene in gene_names:
        url = f"{server}{endpoint}{gene}?"
        response = requests.get(url, headers=headers)
        if response.ok:
            data = response.json()
            gene_ids[gene] = data.get('id')
        else:
            gene_ids[gene] = None

    return gene_ids

if __name__ == "__main__":
    gene_ids = get_ensembl_id(gene_names)
    
    # Create a DataFrame from the gene_ids dictionary
    ensembl_ids_df = pd.DataFrame(list(gene_ids.items()), columns=['Gene Name', 'Ensembl ID'])
    
    # Merge the original DataFrame with the new DataFrame containing Ensembl IDs
    result_df = df.merge(ensembl_ids_df, on='Gene Name', how='left')
    
    # Save the DataFrame to the original CSV file (or you can save to a new file if preferred)
    result_df.to_csv(csv_file_path, index=False)
    
    # Print the Ensembl IDs
    for gene, ensembl_id in gene_ids.items():
        if ensembl_id:
            print(f"The Ensembl ID for {gene} is {ensembl_id}")
        else:
            print(f"Ensembl ID for {gene} not found")
