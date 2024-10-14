from Bio import SeqIO
import numpy as np
import pandas as pd
from tqdm import tqdm

# Function to compute the 20x20 amino acid pair proportion matrix
def compute_adjacent_amino_acid_pair_proportion_matrix(sequences):
    # Initialize a 20x20 zero matrix
    matrix = np.zeros((20, 20), dtype=int)
    total_pairs = 0  # Total number of dipeptides

    for seq in sequences:
        # Traverse each adjacent amino acid pair in the sequence
        for i in range(len(seq) - 1):
            aa1 = seq[i]
            aa2 = seq[i + 1]
            idx1 = aa_index[aa1]
            idx2 = aa_index[aa2]
            # Increment the count for the corresponding dipeptide
            matrix[idx1][idx2] += 1
            total_pairs += 1  # Increment total dipeptide count

    # Convert counts to proportions
    if total_pairs > 0:
        proportion_matrix = matrix / total_pairs
    else:
        proportion_matrix = matrix  # Avoid division by zero

    return proportion_matrix

# Define the 20 standard amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
aa_index = {aa: idx for idx, aa in enumerate(amino_acids)}

# Create a list of all possible dipeptides (400 combinations)
dipeptide_list = [aa1 + aa2 for aa1 in amino_acids for aa2 in amino_acids]

# Initialize a list to collect data
data_list = []

file_list = ['L2_after_i2_cdhit_structured.fasta', 'L2_after_i2_cdhit_disordered.fasta',
             'L3_after_i2_cdhit_structured.fasta', 'L3_after_i2_cdhit_disordered.fasta']

for filename in tqdm(file_list, desc="Processing each file"):
    dataset = './dataset/' + filename
    Level = filename.split('_')[0]
    Structure = filename.split('_')[4].split('.')[0]
    clade_sequences = {}
    # Open and parse the FASTA file
    with open(dataset, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extract the clade information from the header
            header_parts = record.description.strip().split('_')
            if len(header_parts) >= 3:
                clade = header_parts[2]
                # Exclude clades labeled as 'NA'
                if clade != 'NA':
                    # Convert sequence to uppercase and filter non-standard amino acids
                    sequence = ''.join([aa for aa in str(record.seq).upper() if aa in amino_acids])
                    if clade not in clade_sequences:
                        clade_sequences[clade] = []
                    clade_sequences[clade].append(sequence)

    # Compute and collect the proportion data for each clade
    for clade, sequences in clade_sequences.items():
        proportion_matrix = compute_adjacent_amino_acid_pair_proportion_matrix(sequences)

        # Flatten the matrix into a list of proportions in the order of dipeptide_list
        proportions = proportion_matrix.flatten()

        # Create a dictionary with dipeptide proportions
        data = {'Level': Level, 'Structure': Structure, 'Clade': clade}
        for idx, dipeptide in enumerate(dipeptide_list):
            data[dipeptide] = proportions[idx]

        # Append the data to the data_list
        data_list.append(data)

# Create DataFrame from data_list
df = pd.DataFrame(data_list, columns=['Level', 'Structure', 'Clade'] + dipeptide_list)

# Split the DataFrame based on 'Level' and save to separate CSV files
df_L2 = df[df['Level'] == 'L2']
df_L3 = df[df['Level'] == 'L3']

# Save the DataFrames to CSV files
df_L2.to_csv("./L2_dipeptide_proportions.csv", index=False)
df_L3.to_csv("./L3_dipeptide_proportions.csv", index=False)
print("Dipeptide proportion data saved to './L2_dipeptide_proportions.csv' and './L3_dipeptide_proportions.csv'")










###########

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np
import os

# Load the dipeptide proportion data from both files
csv_filename_L2 = "./L2_dipeptide_proportions.csv"
csv_filename_L3 = "./L3_dipeptide_proportions.csv"

df_L2 = pd.read_csv(csv_filename_L2)
df_L3 = pd.read_csv(csv_filename_L3)

# Combine the two DataFrames
df_combined = pd.concat([df_L2, df_L3], ignore_index=True)

# Ensure all columns except 'Level', 'Structure', 'Clade' are numeric
feature_cols = df_combined.columns.difference(['Level', 'Structure', 'Clade'])
df_combined[feature_cols] = df_combined[feature_cols].apply(pd.to_numeric, errors='coerce')

# Handle missing values
df_combined.fillna(0, inplace=True)

# Extract features and labels
X = df_combined[feature_cols].values

# Standardize the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X_scaled)

# Create a DataFrame with the principal components
pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pc_df['Level'] = df_combined['Level']
pc_df['Structure'] = df_combined['Structure']
pc_df['Clade'] = df_combined['Clade']

# Assign colors to 'Level' types
level_color_dict = {'L2': 'blue', 'L3': 'red'}

# Assign markers to 'Clade' types
clades = pc_df['Clade'].unique()
markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', 'h', 'H']  # Extend this list if needed
clade_marker_dict = {clade: markers[i % len(markers)] for i, clade in enumerate(clades)}

# Assign edge colors to 'Structure' types
structure_edgecolor_dict = {'structured': 'black', 'disordered': 'white'}

# Plot the PCA results
plt.figure(figsize=(10, 8))

for idx, row in pc_df.iterrows():
    level = row['Level']
    structure = row['Structure']
    clade = row['Clade']
    x = row['PC1']
    y = row['PC2']
    color = level_color_dict.get(level, 'gray')
    marker = clade_marker_dict.get(clade, 'o')
    edgecolor = structure_edgecolor_dict.get(structure, 'black')
    plt.scatter(x, y, color=color, marker=marker, s=100, edgecolors=edgecolor, linewidths=1)

# Create custom legend entries
import matplotlib.lines as mlines

# Legend entries for 'Level' types (colors)
level_entries = [mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=10, label=level, markerfacecolor=color, markeredgecolor='black') for level, color in level_color_dict.items()]

# Legend entries for 'Clade' types (markers)
clade_entries = [mlines.Line2D([], [], color='gray', marker=clade_marker_dict[clade], linestyle='None', markersize=10, label=clade, markerfacecolor='gray', markeredgecolor='black') for clade in clades]

# Legend entries for 'Structure' types (edge colors)
structure_entries = [mlines.Line2D([], [], color='white', marker='o', linestyle='None', markersize=10, label=structure, markerfacecolor='white', markeredgecolor=edgecolor) for structure, edgecolor in structure_edgecolor_dict.items()]

# Combine legend entries
legend_entries = level_entries + clade_entries + structure_entries

# Add legend to the right-hand side
plt.legend(handles=legend_entries, loc='center left', bbox_to_anchor=(1, 0.5), title='Legend')

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Dipeptide Proportions with Structure Legend (L2 and L3)')
plt.grid(True)
plt.tight_layout()
plt.subplots_adjust(right=0.75)  # Adjust to make room for legend

# Ensure the result directory exists
# os.makedirs('./result/', exist_ok=True)

# Save the plot to a file
plot_filename = "./pca_dipeptide_proportions_L2_L3_structure.png"
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"PCA plot saved to {plot_filename}")
