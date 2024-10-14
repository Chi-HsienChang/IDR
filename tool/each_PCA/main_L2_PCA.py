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

file_list = ['L2_after_i2_cdhit_structured.fasta', 'L2_after_i2_cdhit_disordered.fasta']
# file_list = ['L3_after_i2_cdhit_structured.fasta', 'L3_after_i2_cdhit_disordered.fasta']

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

# Reshape the DataFrame to have dipeptides as rows and clades as columns
# Melt the DataFrame to long format
df_melted = df.melt(id_vars=['Level', 'Structure', 'Clade'], var_name='Dipeptide', value_name='Proportion')

# Create a unique identifier for each clade
df_melted['Clade_ID'] = df_melted['Level'] + '_' + df_melted['Structure'] + '_' + df_melted['Clade']

# Pivot the DataFrame to have dipeptides as rows and Clade_IDs as columns
df_pivot = df_melted.pivot(index='Dipeptide', columns='Clade_ID', values='Proportion')

# Reset index to make 'Dipeptide' a column
df_pivot.reset_index(inplace=True)

# Save the DataFrame to a CSV file
csv_filename = "./result/L2_dipeptide_proportions.csv"
df_pivot.to_csv(csv_filename, index=False)
print(f"Dipeptide proportion data saved to {csv_filename}")

##########

import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np

# Load the dipeptide proportion data
csv_filename = "./result/L2_dipeptide_proportions.csv"  # Ensure this matches your saved file
df = pd.read_csv(csv_filename)

# Set 'Dipeptide' as the index
df.set_index('Dipeptide', inplace=True)

# Transpose the DataFrame so that clades are rows and dipeptides are columns
df_transposed = df.transpose()

# Ensure all column names are strings
df_transposed.columns = df_transposed.columns.astype(str)

# Split the 'Clade_ID' to extract 'Level', 'Structure', and 'Clade'
df_transposed.reset_index(inplace=True)
df_transposed.rename(columns={'index': 'Clade_ID'}, inplace=True)
df_transposed[['Level', 'Structure', 'Clade']] = df_transposed['Clade_ID'].str.split('_', expand=True)

# Reorder columns to move 'Level', 'Structure', 'Clade' to the front
cols = ['Clade_ID', 'Level', 'Structure', 'Clade'] + [col for col in df_transposed.columns if col not in ['Clade_ID', 'Level', 'Structure', 'Clade']]
df_transposed = df_transposed[cols]

# Convert data to numeric, coerce errors to NaN
feature_cols = df_transposed.columns[4:]  # Exclude 'Clade_ID', 'Level', 'Structure', 'Clade'
df_transposed[feature_cols] = df_transposed[feature_cols].apply(pd.to_numeric, errors='coerce')

# Handle missing values (e.g., fill with zeros or drop columns/rows with NaNs)
df_transposed.fillna(0, inplace=True)  # You can choose an appropriate method

# Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(df_transposed[feature_cols])

# Create a DataFrame with the principal components and metadata
pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pc_df = pd.concat([pc_df, df_transposed[['Clade_ID', 'Level', 'Structure', 'Clade']]], axis=1)

# Assign colors to 'Structure' types
structure_types = pc_df['Structure'].unique()
structure_colors = {'structured': 'blue', 'disordered': 'red'}  # Adjust colors if needed

# Assign markers to 'Clade' types
clades = pc_df['Clade'].unique()
markers = [ 'P', '*', 'X', 'h', 'H']  # Extend this list if needed
# markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', 'h', 'H']  # Extend this list if needed
clade_marker_dict = {clade: markers[i % len(markers)] for i, clade in enumerate(clades)}

# Plot the PCA results
plt.figure(figsize=(10, 8))

for idx, row in pc_df.iterrows():
    structure = row['Structure']
    clade = row['Clade']
    x = row['PC1']
    y = row['PC2']
    color = structure_colors.get(structure, 'gray')
    marker = clade_marker_dict.get(clade, 'o')
    plt.scatter(x, y, color=color, marker=marker, s=100)
    
# Create custom legend entries
import matplotlib.lines as mlines

# Legend entries for 'Structure' types (colors)
structure_entries = [mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=10, label=structure) for structure, color in structure_colors.items()]

# Legend entries for 'Clade' types (markers)
clade_entries = [mlines.Line2D([], [], color='black', marker=clade_marker_dict[clade], linestyle='None', markersize=10, label=clade) for clade in clades]

# Combine legend entries
legend_entries = structure_entries + clade_entries

# Add legend to the right-hand side
plt.legend(handles=legend_entries, loc='center left', bbox_to_anchor=(1, 0.5), title='Legend')

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('L2: PCA of Dipeptide Proportions by Structure and Clade')
plt.grid(True)
plt.tight_layout()
plt.subplots_adjust(right=0.75)  # Adjust subplot to make room for the legend

# Save the plot to a file
plot_filename = "./result/L2_pca_dipeptide_proportions.png"
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"PCA plot saved to {plot_filename}")



###########

# import pandas as pd
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt
# import numpy as np

# # Load the dipeptide proportion data
# csv_filename = "./result/dipeptide_proportions.csv"  # Updated filename
# df = pd.read_csv(csv_filename)

# # Set 'Dipeptide' as the index
# df.set_index('Dipeptide', inplace=True)

# # Transpose the DataFrame so that clades are rows and dipeptides are columns
# df_transposed = df.transpose()

# # Ensure all column names are strings
# df_transposed.columns = df_transposed.columns.astype(str)

# # Split the 'Clade_ID' to extract 'Level', 'Structure', and 'Clade'
# df_transposed.reset_index(inplace=True)
# df_transposed.rename(columns={'index': 'Clade_ID'}, inplace=True)
# df_transposed[['Level', 'Structure', 'Clade']] = df_transposed['Clade_ID'].str.split('_', expand=True)

# # Reorder columns to move 'Level', 'Structure', 'Clade' to the front
# cols = ['Clade_ID', 'Level', 'Structure', 'Clade'] + [col for col in df_transposed.columns if col not in ['Clade_ID', 'Level', 'Structure', 'Clade']]
# df_transposed = df_transposed[cols]

# # Convert data to numeric, coerce errors to NaN
# feature_cols = df_transposed.columns[4:]  # Exclude 'Clade_ID', 'Level', 'Structure', 'Clade'
# df_transposed[feature_cols] = df_transposed[feature_cols].apply(pd.to_numeric, errors='coerce')

# # Handle missing values (e.g., fill with zeros or drop columns/rows with NaNs)
# df_transposed.fillna(0, inplace=True)  # You can choose an appropriate method

# # Perform PCA
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(df_transposed[feature_cols])

# # Create a DataFrame with the principal components and metadata
# pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
# pc_df = pd.concat([pc_df, df_transposed[['Clade_ID', 'Level', 'Structure', 'Clade']]], axis=1)

# # Create 'Structure_Clade' combination
# pc_df['Structure_Clade'] = pc_df['Structure'] + '_' + pc_df['Clade']

# # Get unique 'Structure_Clade' combinations
# structure_clade_combinations = pc_df['Structure_Clade'].unique()

# # Assign colors and markers to each 'Structure_Clade' combination
# # Create a colormap with enough colors
# num_combinations = len(structure_clade_combinations)
# colors = plt.cm.get_cmap('tab20', num_combinations).colors
# markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', 'h', 'H', '+', 'x', 'd', '|', '_']  # Extend if needed

# # Ensure we have enough markers
# if num_combinations > len(markers):
#     # Repeat markers if not enough
#     markers = markers * (num_combinations // len(markers) + 1)

# structure_clade_style_dict = {}
# for i, combination in enumerate(structure_clade_combinations):
#     color = colors[i % len(colors)]
#     marker = markers[i]
#     structure_clade_style_dict[combination] = {'color': color, 'marker': marker}

# # Plot the PCA results
# plt.figure(figsize=(10, 8))

# for idx, row in pc_df.iterrows():
#     combination = row['Structure_Clade']
#     x = row['PC1']
#     y = row['PC2']
#     style = structure_clade_style_dict.get(combination, {'color': 'black', 'marker': 'o'})
#     plt.scatter(x, y, color=style['color'], marker=style['marker'], s=100, label=combination)

# # Remove duplicate labels in legend
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))

# # Add legend to the right-hand side
# plt.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5), title='Structure_Clade')

# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.title('PCA of Dipeptide Proportions by Structure and Clade')
# plt.grid(True)
# plt.tight_layout()
# plt.subplots_adjust(right=0.75)  # Adjust subplot to make room for the legend

# # Save the plot to a file
# plot_filename = "./result/pca_dipeptide_proportions.png"
# plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
# plt.show()
# print(f"PCA plot saved to {plot_filename}")



###########

# import pandas as pd
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt
# import numpy as np

# # Load the dipeptide proportion data
# csv_filename = "./result/dipeptide_proportions.csv"  # Updated filename
# df = pd.read_csv(csv_filename)

# # Set 'Dipeptide' as the index
# df.set_index('Dipeptide', inplace=True)

# # Transpose the DataFrame so that clades are rows and dipeptides are columns
# df_transposed = df.transpose()

# # Ensure all column names are strings
# df_transposed.columns = df_transposed.columns.astype(str)

# # Convert the DataFrame to numeric, coerce errors to NaN
# df_transposed = df_transposed.apply(pd.to_numeric, errors='coerce')

# # Handle missing values (e.g., fill with zeros or drop columns/rows with NaNs)
# df_transposed.fillna(0, inplace=True)  # You can choose an appropriate method

# # Perform PCA
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(df_transposed)

# # Create a DataFrame with the principal components
# pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
# pc_df['Clade_ID'] = df_transposed.index

# # Plot the PCA results
# plt.figure(figsize=(8, 6))
# plt.scatter(pc_df['PC1'], pc_df['PC2'], color='blue')

# # Annotate each point with the clade ID
# for i, clade_id in enumerate(pc_df['Clade_ID']):
#     plt.annotate(clade_id, (pc_df['PC1'][i], pc_df['PC2'][i]))

# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.title('PCA of Dipeptide Proportions by Clade')
# plt.grid(True)
# plt.tight_layout()

# # Save the plot to a file
# plot_filename = "./result/pca_dipeptide_proportions.png"
# plt.savefig(plot_filename, dpi=300)
# plt.show()
# print(f"PCA plot saved to {plot_filename}")

