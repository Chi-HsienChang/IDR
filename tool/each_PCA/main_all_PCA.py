import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# File paths
file1 = './L2_dipeptide_proportions.csv'
file2 = './L3_dipeptide_proportions.csv'

# Load the CSV files
df_L2 = pd.read_csv(file1)
df_L3 = pd.read_csv(file2)

# Get observation labels
df_L2_obs = df_L2.columns[1:].tolist()  # Exclude 'Dipeptide' column
df_L3_obs = df_L3.columns[1:].tolist()

# Combine the two dataframes
combined_df = pd.concat([df_L2, df_L3], axis=0).reset_index(drop=True)

# Handle missing values: Fill NaNs with 0
combined_df = combined_df.fillna(0)

# Transpose the data
transposed_df = combined_df.set_index('Dipeptide').T

# Ensure all column names are strings
transposed_df.columns = transposed_df.columns.astype(str)

# Map observations to groups
obs_to_group = {obs: 'L2' for obs in df_L2_obs}
obs_to_group.update({obs: 'L3' for obs in df_L3_obs})

# Perform PCA on the transposed data
pca = PCA(n_components=2)
principal_components = pca.fit_transform(transposed_df)

# Get variance ratios
variance_ratios = pca.explained_variance_ratio_ * 100  # Convert to percentage

# Create a DataFrame for the PCA results
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pca_df['Observation'] = transposed_df.index
pca_df['Group'] = pca_df['Observation'].map(obs_to_group)

# Save PCA results to a CSV file
pca_df.to_csv('./PCA_transposed_results.csv', index=False)

# Visualize the PCA results
plt.figure(figsize=(10, 7))

# Define colors
colors = {'L2': 'blue', 'L3': 'red'}

# Increase font size
plt.rcParams.update({'font.size': 14})

# Plot each group
for group, color in colors.items():
    idx = pca_df['Group'] == group
    plt.scatter(pca_df.loc[idx, 'PC1'], pca_df.loc[idx, 'PC2'], alpha=0.7, label=group, color=color)
    for i in pca_df.loc[idx].index:
        plt.annotate(pca_df['Observation'][i], (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=12, alpha=0.7)

# Add variance explained to axis labels
# plt.title('PCA of Observations (Transposed Data)')
plt.xlabel(f'Principal Component 1 ({variance_ratios[0]:.2f}% Variance)')
plt.ylabel(f'Principal Component 2 ({variance_ratios[1]:.2f}% Variance)')
plt.grid(True)
plt.legend()

print("PCA analysis completed and results saved to 'PCA_transposed_results.csv'.")
plt.savefig("PCA_1118_3.png")






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# import pandas as pd
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt

# # File paths
# file1 = './L2_dipeptide_proportions.csv'
# file2 = './L3_dipeptide_proportions.csv'

# # Load the CSV files
# df_L2 = pd.read_csv(file1)
# df_L3 = pd.read_csv(file2)

# # Get observation labels
# df_L2_obs = df_L2.columns[1:].tolist()  # Exclude 'Dipeptide' column
# df_L3_obs = df_L3.columns[1:].tolist()

# # Combine the two dataframes
# combined_df = pd.concat([df_L2, df_L3], axis=0).reset_index(drop=True)

# # Handle missing values: Fill NaNs with 0
# combined_df = combined_df.fillna(0)

# # Transpose the data
# transposed_df = combined_df.set_index('Dipeptide').T

# # Ensure all column names are strings
# transposed_df.columns = transposed_df.columns.astype(str)

# # Map observations to groups
# obs_to_group = {obs: 'L2' for obs in df_L2_obs}
# obs_to_group.update({obs: 'L3' for obs in df_L3_obs})

# # Perform PCA on the transposed data
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(transposed_df)

# # Create a DataFrame for the PCA results
# pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
# pca_df['Observation'] = transposed_df.index
# pca_df['Group'] = pca_df['Observation'].map(obs_to_group)

# # Save PCA results to a CSV file
# pca_df.to_csv('./PCA_transposed_results.csv', index=False)

# # Visualize the PCA results
# plt.figure(figsize=(10, 7))

# # Define colors
# colors = {'L2': 'blue', 'L3': 'red'}

# # Increase font size
# plt.rcParams.update({'font.size': 14})

# # Plot each group
# for group, color in colors.items():
#     idx = pca_df['Group'] == group
#     plt.scatter(pca_df.loc[idx, 'PC1'], pca_df.loc[idx, 'PC2'], alpha=0.7, label=group, color=color)
#     for i in pca_df.loc[idx].index:
#         plt.annotate(pca_df['Observation'][i], (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=12, alpha=0.7)

# # plt.title('PCA of Observations (Transposed Data)')
# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.grid(True)
# plt.legend()

# print("PCA analysis completed and results saved to 'PCA_transposed_results.csv'.")
# plt.savefig("PCA_1118_2.png")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# import pandas as pd
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt

# # File paths
# file1 = './L2_dipeptide_proportions.csv'
# file2 = './L3_dipeptide_proportions.csv'

# # Load the CSV files
# df_L2 = pd.read_csv(file1)
# df_L3 = pd.read_csv(file2)

# # Combine the two dataframes
# combined_df = pd.concat([df_L2, df_L3], axis=0).reset_index(drop=True)

# # Handle missing values: Fill NaNs with 0 (or another strategy like mean/median)
# combined_df = combined_df.fillna(0)

# # Transpose the data
# transposed_df = combined_df.set_index('Dipeptide').T

# # Ensure all column names are strings
# transposed_df.columns = transposed_df.columns.astype(str)

# # Perform PCA on the transposed data
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(transposed_df)

# # Create a DataFrame for the PCA results
# pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
# pca_df['Observation'] = transposed_df.index

# # Save PCA results to a CSV file
# pca_df.to_csv('./PCA_transposed_results.csv', index=False)

# # Visualize the PCA results
# plt.figure(figsize=(10, 7))
# plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
# for i, label in enumerate(pca_df['Observation']):
#     plt.annotate(label, (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=8, alpha=0.7)
# plt.title('PCA of Observations (Transposed Data)')
# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.grid(True)

# print("PCA analysis completed and results saved to 'PCA_transposed_results.csv'.")
# plt.savefig("PCA_1118.png")
