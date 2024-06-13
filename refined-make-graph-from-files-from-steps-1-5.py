import pandas as pd
import matplotlib.pyplot as plt
import re

# Load the exon data
exon_data = pd.read_csv('exon_data.csv')

# Load the new mutation data
new_mutation_data = pd.read_csv('mutation_data.csv')

# Extract chromosome and position information from the Genomic Location column
def parse_genomic_location(location):
    if isinstance(location, str):
        match = re.match(r'NC_(\d{6}).(\d+):g.(\d+)[A-Z]>', location)
        if match:
            return int(match.group(1)), int(match.group(3))
    return None, None

# Apply the parsing function to the new mutation_data dataframe
new_mutation_data[['Chromosome', 'Position']] = new_mutation_data['Genomic Location'].apply(
    lambda x: pd.Series(parse_genomic_location(x))
)

# Remove rows with parsing errors
new_mutation_data = new_mutation_data.dropna(subset=['Chromosome', 'Position'])

# Convert the chromosome data to string to match the exon data format
new_mutation_data['Chromosome'] = new_mutation_data['Chromosome'].astype(int).astype(str)

# Map mutations to exons
def map_mutation_to_exon(row, exon_data):
    exon = exon_data[(exon_data['Seq Region Name'] == row['Chromosome']) &
                     (exon_data['Start'] <= row['Position']) &
                     (exon_data['End'] >= row['Position'])]
    if not exon.empty:
        return exon.iloc[0]['Exon ID']
    return None

new_mutation_data['Exon ID'] = new_mutation_data.apply(lambda row: map_mutation_to_exon(row, exon_data), axis=1)

# Filter out mutations that do not map to any exon
mapped_mutations = new_mutation_data.dropna(subset=['Exon ID'])

# Plot the distribution of mutations across exons
mutation_counts = mapped_mutations['Exon ID'].value_counts().sort_index()

plt.figure(figsize=(12, 6))
mutation_counts.plot(kind='bar')
plt.title('Distribution of Mutations Across Exons')
plt.xlabel('Exon ID')
plt.ylabel('Number of Mutations')
plt.xticks(rotation=90)
plt.show()
