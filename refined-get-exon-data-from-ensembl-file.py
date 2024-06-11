import json

# Load the JSON file
file_path = 'AIRE_ensembl.json'

# Load the JSON data from the file
with open(file_path, 'r') as file:
    data = json.load(file)

# Extract and print the start, end, and seq_region_name under the Exon key in Transcript
for transcript in data.get('Transcript', []):
    if 'Exon' in transcript:
        print(f"Transcript ID: {transcript.get('id', 'N/A')}")
        for exon in transcript['Exon']:
            start = exon.get('start', 'N/A')
            end = exon.get('end', 'N/A')
            seq_region_name = exon.get('seq_region_name', 'N/A')
            print(f"  Exon ID: {exon.get('id', 'N/A')} - Start: {start}, End: {end}, Seq Region Name: {seq_region_name}")
        print()  # Add a blank line between transcripts
