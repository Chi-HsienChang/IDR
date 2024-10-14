from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

# Initialize a dictionary to hold sequences by clade


file_list = ['L2_after_i2_cdhit_structured.fasta', 'L3_after_i2_cdhit_structured.fasta', 'L2_after_i2_cdhit_disordered.fasta', 'L3_after_i2_cdhit_disordered.fasta']

for file in tqdm(file_list, "each file"):
    dataset = './dataset/'+file
    clade_sequences = {}
    # Open and parse the FASTA file
    with open(dataset, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extract the clade information from the header
            header_parts = record.description.strip().split('_')
            if len(header_parts) >= 3:
                clade = header_parts[2]
                # Convert sequence to uppercase and filter non-standard amino acids
                sequence = ''.join([aa for aa in str(record.seq).upper() if aa in amino_acids])
                if clade not in clade_sequences:
                    clade_sequences[clade] = []
                clade_sequences[clade].append(sequence)

        #############
        # Compute and save the sorted proportion lists and plots for each clade
        for clade, sequences in clade_sequences.items():
            proportion_matrix = compute_adjacent_amino_acid_pair_proportion_matrix(sequences)

            # Flatten the matrix into a list of tuples: (dipeptide, proportion)
            dipeptides = []
            for i in range(20):
                for j in range(20):
                    dipeptide = amino_acids[i] + amino_acids[j]
                    proportion = proportion_matrix[i][j]
                    dipeptides.append((dipeptide, proportion))

            # Sort the list from highest to lowest proportion
            dipeptides_sorted = sorted(dipeptides, key=lambda x: x[1], reverse=True)

            # Convert to DataFrame
            df = pd.DataFrame(dipeptides_sorted, columns=['Dipeptide', 'Proportion'])

            # Save to CSV file
            # Optionally, you can translate clade names if needed
            clade_translation = {
                'Stramenopiles': 'Stramenopiles',
                'Metazoa': 'Metazoa',
                'Fungi': 'Fungi',
                # Add more translations if needed
            }
            clade_name = clade_translation.get(clade, clade)
            csv_filename = f"./result/{file.split('_')[0]}/csv/{file.split('_')[0]}_{file.split('_')[4].split('.')[0]}_{clade_name}.csv"
            df.to_csv(csv_filename, index=False)
            print(f"Proportion list for clade {clade} saved to {csv_filename}")

            # Plot the top 20 dipeptides as a bar chart
            top_n = 20
            top_dipeptides = dipeptides_sorted[:top_n]
            dipeptide_names = [dp[0] for dp in top_dipeptides]
            proportions = [dp[1] for dp in top_dipeptides]

            plt.figure(figsize=(12, 6))
            plt.bar(dipeptide_names, proportions, color='skyblue')
            plt.xlabel('Dipeptide')
            plt.ylabel('Proportion')
            plt.title(f'Top {top_n} Dipeptides in Clade {clade_name}')
            plt.xticks(rotation=90)
            plt.tight_layout()

            # Save plot to file
            plot_filename = f"./result/{file.split('_')[0]}/png/{file.split('_')[0]}_{file.split('_')[4].split('.')[0]}_{clade_name}.png"
            plt.savefig(plot_filename, dpi=300)
            plt.close()  # Close the figure to free memory
            print(f"Bar chart for clade {clade} saved to {plot_filename}")
