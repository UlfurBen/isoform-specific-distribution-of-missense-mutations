import json

# Define the file path
file_path = 'O43918_refined.json'

# Load the JSON data from the file
with open(file_path, 'r') as file:
    data = json.load(file)

# Extract the sequence value under the sequence key and print them
for i, entry in enumerate(data, start=1):
    if 'sequence' in entry:
        print(f"Sequence {i}:")
        print(entry['sequence']['sequence'])
        print()

