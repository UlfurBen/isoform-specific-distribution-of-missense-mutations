import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

# Sample mutation data
data = {
    'Mutation': ['c.16442G>A', 'c.16443A>G', 'c.12000T>C', 'c.12500G>T'],
    'Protein_Change': ['p.Cys5481Tyr', 'p.Arg5482Gly', 'p.Leu4000Pro', 'p.Gly4200Val'],
    'Exon': [48, 48, 32, 34]  # Example exon numbers
}

# Create DataFrame
df = pd.DataFrame(data)

# Count mutations by exon
mutation_counts = df['Exon'].value_counts().sort_index()

# Plot mutation counts
plt.bar(mutation_counts.index, mutation_counts.values)
plt.xlabel('Exon Number')
plt.ylabel('Mutation Count')
plt.title('Mutation Enrichment Across Exons of KMT2D')
plt.show()

# Statistical test (example)
observed = mutation_counts.values
expected = [len(df) / len(mutation_counts)] * len(mutation_counts)  # assuming uniform distribution

chi2, p_value = chi2_contingency([observed, expected])[:2]
print(f'Chi-square statistic: {chi2}, p-value: {p_value}')
