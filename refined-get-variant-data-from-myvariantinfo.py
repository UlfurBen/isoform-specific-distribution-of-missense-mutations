import requests
import pandas as pd

# List of genes to fetch mutation information for
genes = ["BRCA1", "BRCA2", "TP53"]  # Replace with your list of genes

# Function to get mutation information for a single gene using MyVariant.info API
def get_gene_mutation_info(gene):
    url = f"https://myvariant.info/v1/query?q={gene}&fields=all"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return None

# Function to filter missense mutations from the results
def filter_missense_mutations(mutations):
    missense_mutations = []
    for mutation in mutations.get('hits', []):
        if 'snpeff' in mutation:
            for annotation in mutation['snpeff'].get('ann', []):
                if isinstance(annotation, dict) and annotation.get('effect') == 'missense_variant':
                    missense_mutations.append({
                        "gene": mutation.get("query", ""),
                        "mutation_id": mutation.get("_id", ""),
                        "effect": annotation.get('effect', ''),
                        "impact": annotation.get('putative_impact', ''),
                        "hgvs_c": annotation.get('hgvs_c', ''),
                        "hgvs_p": annotation.get('hgvs_p', ''),
                        "variant_allele": annotation.get('variant_allele', ''),
                        "transcript_id": annotation.get('transcript_id', ''),
                        "gene_id": annotation.get('gene_id', ''),
                    })
    return missense_mutations

# Collect missense mutation information for all genes
gene_mutation_data = []

for gene in genes:
    result = get_gene_mutation_info(gene)
    if result:
        filtered_data = filter_missense_mutations(result)
        gene_mutation_data.extend(filtered_data)

# Convert the data into a pandas DataFrame
df = pd.DataFrame(gene_mutation_data)

# Save the DataFrame to a CSV file
output_csv_file = "missense_gene_mutation_data.csv"
df.to_csv(output_csv_file, index=False)

print(f"Missense mutation data saved to {output_csv_file}")
