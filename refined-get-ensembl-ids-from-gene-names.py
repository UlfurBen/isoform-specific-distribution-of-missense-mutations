import os
import requests
import json
import pandas as pd

csv_file_path = "/Users/ulfurfjolnisson/Downloads/The Epigenetic Machinery.csv"

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
    result_df = pd.DataFrame(list(gene_ids.items()), columns=['Gene Name', 'Ensembl ID'])
    
    # Save the DataFrame to a CSV file
    result_df.to_csv('ensembl_gene_ids.csv', index=False)
    
    for gene, ensembl_id in gene_ids.items():
        if ensembl_id:
            print(f"The Ensembl ID for {gene} is {ensembl_id}")
        else:
            print(f"Ensembl ID for {gene} not found")
