import json
import csv

# Load the JSON file
file_path = 'AIRE_ensembl.json'

# Load the JSON data from the file
with open(file_path, 'r') as file:
    data = json.load(file)

# Open a CSV file to write the data
with open('transcripts_exons.csv', 'w', newline='') as csvfile:
    fieldnames = ['Transcript ID', 'Exon ID', 'Start', 'End', 'Seq Region Name']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()

    # Extract and write the start, end, and seq_region_name under the Exon key in Transcript
    for transcript in data.get('Transcript', []):
        transcript_id = transcript.get('id', 'N/A')
        if 'Exon' in transcript:
            for exon in transcript['Exon']:
                start = exon.get('start', 'N/A')
                end = exon.get('end', 'N/A')
                seq_region_name = exon.get('seq_region_name', 'N/A')
                exon_id = exon.get('id', 'N/A')

                writer.writerow({
                    'Transcript ID': transcript_id,
                    'Exon ID': exon_id,
                    'Start': start,
                    'End': end,
                    'Seq Region Name': seq_region_name
                })
