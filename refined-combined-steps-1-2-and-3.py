import os
import requests
import json
import pandas as pd

csv_file_path = "/Users/ulfurfjolnisson/Work/The Epigenetic Machinery.csv"

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
    result_df = df.merge(ensembl_ids_df, on='Gene Name', how='left', suffixes=('', '_new'))
    
    # Save the DataFrame to the original CSV file (or you can save to a new file if preferred)
    result_df.to_csv(csv_file_path, index=False)
    
    # Print the Ensembl IDs
    for gene, ensembl_id in gene_ids.items():
        if ensembl_id:
            print(f"The Ensembl ID for {gene} is {ensembl_id}")
        else:
            print(f"Ensembl ID for {gene} not found")

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






# Directory containing JSON files
directory = 'isoform_information'

# Initialize a list to store the extracted data
exon_data = []

# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.json'):
        file_path = os.path.join(directory, filename)
        
        # Load the JSON data from the file
        with open(file_path, 'r') as file:
            data = json.load(file)
        
        # Extract the start, end, and seq_region_name under the Exon key in Transcript
        for transcript in data.get('Transcript', []):
            transcript_id = transcript.get('id', 'N/A')
            if 'Exon' in transcript:
                for exon in transcript['Exon']:
                    chr = exon.get('seq_region_name', 'N/A')
                    exon_id = exon.get('id', 'N/A')
                    start = exon.get('start', 'N/A')
                    end = exon.get('end', 'N/A')
                    seq_region_name = exon.get('seq_region_name', 'N/A')
                    exon_data.append([transcript_id, exon_id, start, end, seq_region_name])

# Create a DataFrame from the list
exon_df = pd.DataFrame(exon_data, columns=['Transcript ID', 'Exon ID', 'Start', 'End', 'Seq Region Name'])

# Save the DataFrame to a CSV file
exon_df.to_csv('exon_data.csv', index=False)

print("Data saved to exon_data.csv")
