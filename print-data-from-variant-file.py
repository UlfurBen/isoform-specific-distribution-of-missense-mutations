import json

# Load the JSON file
file_path = 'Q9NP73_variation.json'

# Load the JSON data from the file
with open(file_path, 'r') as file:
    data = json.load(file)

# Iterate through the JSON data to extract and print the required information
for feature in data.get('features', []):
    consequence_type = feature.get('consequenceType', 'N/A')
    genomic_location = ', '.join(feature.get('genomicLocation', []))
    
    population_frequencies = feature.get('populationFrequencies', [])
    if population_frequencies:
        for frequency_data in population_frequencies:
            frequency = frequency_data.get('frequency', 'N/A')
            print(f"Consequence Type: {consequence_type}, Genomic Location: {genomic_location}, Frequency: {frequency}")
    else:
        print(f"Consequence Type: {consequence_type}, Genomic Location: {genomic_location}, Frequency: N/A")
