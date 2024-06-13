import json
import csv
import os

# Define the CSV file path
csv_file_path = 'accession_ids.csv'

# Directory containing the JSON files
directory = 'isoforms'

# Open the CSV file for writing
with open(csv_file_path, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csvwriter = csv.writer(csvfile)
    
    # Write the header
    csvwriter.writerow(['Accession ID', 'Accession', 'Sequence'])

    for filename in os.listdir(directory):
        if filename.endswith('.json'):
            file_path = os.path.join(directory, filename)
            
            # Load the JSON data from the file
            with open(file_path, 'r') as file:
                data = json.load(file)
            
            # Check if the data is a list or dictionary and process accordingly
            if isinstance(data, list):
                for entry in data:
                    if isinstance(entry, dict):
                        accession = entry.get('accession', 'N/A')
                        sequence = entry.get('sequence', {}).get('sequence', 'N/A')
                        csvwriter.writerow([filename, accession, sequence])
            elif isinstance(data, dict):
                for feature in data.get('features', []):
                    accession = feature.get('accession', 'N/A')
                    sequence = feature.get('sequence', {}).get('sequence', 'N/A')
                    csvwriter.writerow([filename, accession, sequence])
            else:
                print(f"Unexpected data format in file: {file_path}")
