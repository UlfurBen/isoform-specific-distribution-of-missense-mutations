import requests
import json

# URL of the UniProt entry in JSON format for KMT2D
id = 'A0A8I5QJJ1'
url = f"https://www.ebi.ac.uk/proteins/api/variation/{id}?format=json"




# Send a request to the UniProt API
response = requests.get(url)

data = response.json()
print(data)
# Check if the request was successful

