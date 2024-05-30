import requests
import json

uniprot_id = "O43918"

# URL of the UniProt entry in JSON format
url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"

# Send a request to the UniProt API
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    data = response.json()
    
    # Initialize a list to hold isoform and variant information
    isoform_data = []
    
    # Extract isoform variant information
    if 'comments' in data:
        for comment in data['comments']:
            if comment['commentType'] == 'ALTERNATIVE PRODUCTS':
                for isoform in comment['isoforms']:
                    isoform_entry = {
                        "isoform_id": isoform['isoformIds'],
                        "name": isoform['name']['value'],
                        "sequence_ids": isoform.get('sequenceIds', [])
                    }
                    
                    # Initialize a list to hold variant pathogenicity information
                    variants_info = []
                    
                    # Extract variant information from features
                    for feature in data.get('features', []):
                        if feature['type'] == 'variant':
                            variant_info = {
                                "location": feature['location']['start'],
                                "description": feature.get('description', ''),
                                "pathogenicity": feature.get('evidences', [])
                            }
                            variants_info.append(variant_info)
                    
                    # Add variants information to the isoform entry
                    isoform_entry["variants"] = variants_info
                    
                    # Add the isoform entry to the main list
                    isoform_data.append(isoform_entry)
    
    # Save the combined information to a JSON file
    with open('isoforms.json', 'w') as f:
        json.dump(isoform_data, f, indent=4)
    
    print("Isoform variant information has been saved to isoforms.json")
else:
    print(f"Failed to retrieve data from UniProt. Status code: {response.status_code}")
