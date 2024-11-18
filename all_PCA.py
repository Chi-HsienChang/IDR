from Bio import SeqIO
import numpy as np
import pandas as pd
from tqdm import tqdm

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

# Calculate explained variance
explained_variance = pca.explained_variance_ratio_

# Create a DataFrame with the principal components
pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pc_df['Level'] = df_combined['Level']
pc_df['Structure'] = df_combined['Structure']
pc_df['Clade'] = df_combined['Clade']

# Assign colors to 'Level' types
level_color_dict = {'L2': 'blue', 'L3': 'red'}

# Assign markers to 'Clade' types
clades = pc_df['Clade'].unique()
markers = ['s', '^', 'D', 'v', 'P', '*', 'X', 'h', 'H']  # Extend this list if needed
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

# Add explained variance information
plt.text(1.05, 1.02, f'Explained Variance:\nPC1: {explained_variance[0]:.2%}\nPC2: {explained_variance[1]:.2%}', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')

plt.tight_layout()
plt.subplots_adjust(right=0.75)  # Adjust to make room for legend

# Ensure the result directory exists
# os.makedirs('./result/', exist_ok=True)

# Save the plot to a file
plot_filename = "./pca_dipeptide_proportions_L2_L3_structure_IDR_1118.png"
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
# plt.show()
print(f"PCA plot saved to {plot_filename}")
