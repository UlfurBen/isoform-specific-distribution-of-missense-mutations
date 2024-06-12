import json
import csv

# Define the file paths
json_file_path = 'O43918_isoform_refined.json'
csv_file_path = 'accession_ids.csv'

# Load the JSON data from the file
with open(json_file_path, 'r') as file:
    data = json.load(file)

# Open the CSV file for writing
with open(csv_file_path, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csvwriter = csv.writer(csvfile)
    
    # Write the header
    csvwriter.writerow(['accession id', 'accession'])

    # Extract the sequence value under the sequence key and write them to the CSV file
    for i, entry in enumerate(data, start=1):
        if 'accession' in entry:
            accession = entry['accession']
            csvwriter.writerow([f"Accession {i}", accession])
