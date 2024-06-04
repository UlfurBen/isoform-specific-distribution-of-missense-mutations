import requests
import sys
import json

# Define the UniProt API URL
uniprot_api = 'https://rest.uniprot.org/uniprotkb/'

def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response

def get_isoform_info(gene_id):
    # Construct the API request for isoform information
    api_request = f'{uniprot_api}{gene_id}.json'

    # Send the request
    response = get_url(api_request)

    # Parse the JSON response
    data = response.json()

    # Print the API response
    print("API Response:")
    print(json.dumps(data, indent=4))

    # Check if the response contains isoform information
    if 'features' in data:
        isoforms = [feature for feature in data['features'] if feature['type'] == 'CHAIN']

        # Print isoform information
        for isoform in isoforms:
            print(f"Isoform: {isoform['id']}")
            print(f"Sequence: {isoform['sequence']}")
            print(f"Length: {len(isoform['sequence'])}")
            print()
    else:
        print("No isoform information found.")

# Example usage
gene_id = 'O14686'
get_isoform_info(gene_id)
