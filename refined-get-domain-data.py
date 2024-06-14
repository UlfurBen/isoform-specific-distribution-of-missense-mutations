import csv
import requests
from bs4 import BeautifulSoup
import pandas as pd

# Read the input CSV file into a DataFrame
input_csv_path = "The Epigenetic Machinery.csv"
df = pd.read_csv(input_csv_path)

# Extract the Uniprot IDs into a list
uniprot_ids = df['UniProt'].tolist()

# Loop over each Uniprot ID and process the data
for uniprot_id in uniprot_ids:
    # URL of the web page with the current Uniprot ID
    url = f"https://web.expasy.org/cgi-bin/protparam/protparam?{uniprot_id}"

    # Send a request to fetch the web page
    response = requests.get(url)
    if response.status_code == 200:
        # Parse the web page content
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Extracting all data within <pre> tags
        pre_tags = soup.find_all('pre')
        
        # Save the extracted data to a file
        output_file_path = f'domain_info_{uniprot_id}.txt'
        with open(output_file_path, 'w') as file:
            for pre in pre_tags:
                file.write(pre.text + '\n')
        
        # Read data from the text file
        with open(output_file_path, 'r') as file:
            lines = file.readlines()

        # Process the lines and extract relevant columns
        data = []
        for line in lines:
            if line.startswith('FT'):
                parts = line.split()
                feature_type = parts[1]
                range_ = parts[2]
                description = ' '.join(parts[3:])
                
                # Add a single quote to prevent misinterpretation of numbers as dates
                range_ = f"'{range_}"
                
                data.append([feature_type, range_, description])

        # Write the processed data to a CSV file
        csv_output_path = f'domain_info_{uniprot_id}.csv'
        with open(csv_output_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Feature Type', 'Range', 'Description'])  # Write header
            csvwriter.writerows(data)  # Write data rows

        print(f"Data for {uniprot_id} has been written to {csv_output_path}")
    else:
        print(f"Failed to retrieve the page for {uniprot_id}. Status code: {response.status_code}")
