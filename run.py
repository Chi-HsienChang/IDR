from Bio import SeqIO
import numpy as np
import pandas as pd


amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
aa_index = {aa: idx for idx, aa in enumerate(amino_acids)}


clade_sequences = {}

dataset = './dataset/L2_after_i2_cdhit.fasta'

with open(dataset, 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
    
        header_parts = record.description.strip().split('_')
        if len(header_parts) >= 3:
            clade = header_parts[2]
        
            sequence = ''.join([aa for aa in str(record.seq).upper() if aa in amino_acids])
            if clade not in clade_sequences:
                clade_sequences[clade] = []
            clade_sequences[clade].append(sequence)


def compute_adjacent_amino_acid_pair_matrix(sequences):
  
    matrix = np.zeros((20, 20), dtype=int)
    
    for seq in sequences:
       
        for i in range(len(seq) - 1):
            aa1 = seq[i]
            aa2 = seq[i + 1]
            idx1 = aa_index[aa1]
            idx2 = aa_index[aa2]
            
            matrix[idx1][idx2] += 1
    return matrix


for clade, sequences in clade_sequences.items():
    matrix = compute_adjacent_amino_acid_pair_matrix(sequences)
    df = pd.DataFrame(matrix, index=amino_acids, columns=amino_acids)
    print(f"Amino acid pair count matrix for clade {clade}:")
    print(df)
    print("\n")
