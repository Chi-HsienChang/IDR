import pandas as pd

# File paths
fasta_file_path = './aligned_real.fasta'
tsv_file_path = './LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_lineage = {}
name_to_taxID = {}


# Process the FASTA file
with open(fasta_file_path, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            # Extract the identifier and remove the leading '>'
            identifier = line[1:]  # Remove '>' from the start of the header
            parts = identifier.split('|')
            name = parts[0]  # Assume the first part of the header is the name, now cleaned of '>'
            taxID = parts[-2].split(':')[-1]  # Extract taxID from the last part

            # Retrieve the Named Lineage using taxID
            named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
            # Map the name to its corresponding Named Lineage
            name_to_lineage[name] = named_lineage
            name_to_taxID[name] = taxID




from Bio import SeqIO

# File paths
input_fasta = 'L2_after_i5_cdhit_IDR.fasta'
output_fasta = input_fasta.split('.')[0] + '_modified.fasta'

# # Example mappings
# name_to_taxID = {'seq1': '12345', 'seq2': '67890'}
# name_to_lineage = {'seq1': 'CladeA', 'seq2': 'CladeB'}

# Processing the FASTA file
with open(output_fasta, 'w') as outfile:
    for record in SeqIO.parse(input_fasta, 'fasta'):
        original_name = record.id.split('|')[0]
        taxID = name_to_taxID.get(original_name, 'unknown_taxID')
        clade = name_to_lineage.get(original_name, 'unknown_clade')

        # Update header
        new_name = f"{original_name}_{input_fasta.split('.')[0].split('_')[0]}_{clade.split(',')[4]}_{taxID}"
        record.id = new_name
        record.description = ""  # Clear description to avoid duplication
        SeqIO.write(record, outfile, 'fasta')

print(f"Modified FASTA file saved to {output_fasta}")
