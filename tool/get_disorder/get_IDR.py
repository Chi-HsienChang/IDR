from Bio import SeqIO
from ipdb import set_trace


# Load all.fasta into a dictionary mapping IDs to records (including headers)
all_sequences = {}
with open('../../dataset/interpro_structured_disordered.fasta', 'r') as all_file:
    for record in SeqIO.parse(all_file, 'fasta'):
        # Extract ID up to the first '|' character
        id_all = record.id.split('|')[0]
        all_sequences[id_all] = record  # Store the entire record

# Load target.fasta into a list of tuples (ID, sequence, header)
target_sequences = []
file = 'L2_type_iteration_5_after_cluster.fasta'
type = file.split('_')[0]
with open('../../dataset/'+file, 'r') as target_file:
    for record in SeqIO.parse(target_file, 'fasta'):
        # Extract ID up to the first '_' character
        id_target = record.id.split('_')[0]
        
        # Keep the original header from target.fasta
        header_target = record.description
        # Convert the sequence to uppercase
        seq_target = str(record.seq).upper()
        target_sequences.append((id_target, seq_target, header_target))


# set_trace()
# For each target sequence, find the corresponding sequence in all.fasta
# and extract the sequence to the right of the target sequence
results = []
for id_target, seq_target, header_target in target_sequences:
    if id_target in all_sequences:
        record_all = all_sequences[id_target]
        # set_trace()
        seq_all = str(record_all.seq).upper()
        # Find where seq_target occurs in seq_all
        index = seq_all.find(seq_target)
        if index != -1:
            # Extract sequence to the right of the target sequence
            start = index + len(seq_target)
            right_seq = seq_all[start:]
            # Use the original header from target.fasta
            header = header_target
            results.append((header, right_seq))
        else:
            print(f"Target sequence not found in all.fasta sequence for ID {id_target}")
            # Handle this case (store empty sequence with target header)
            header = header_target
            results.append((header, ''))
    else:
        print(f"ID {id_target} not found in all.fasta")
        # Handle this case (store empty sequence with target header)
        header = header_target
        results.append((header, ''))

# Write the results to a new FASTA file, keeping the original target headers
with open(f'{type}_after_i5_cdhit_IDR.fasta', 'w') as output_file:
    for header, seq_result in results:
        output_file.write(f">{header}\n")  # Use the header from target.fasta
        output_file.write(f"{seq_result}\n")

####

# from Bio import SeqIO

# # Load all.fasta into a dictionary mapping IDs to sequences
# all_sequences = {}
# with open('../dataset/test_structured_disordered.fasta', 'r') as all_file:
#     for record in SeqIO.parse(all_file, 'fasta'):
#         # Extract ID up to the first '|' character
#         id_all = record.id.split('|')[0]
#         all_sequences[id_all] = str(record.seq)
#         print(all_sequences)

# # Load target.fasta into a list of tuples (ID, sequence)
# target_sequences = []
# with open('../dataset/test_L2_after_i2_cdhit.fasta', 'r') as target_file:
#     for record in SeqIO.parse(target_file, 'fasta'):
#         # Extract ID up to the first '_' character
#         id_target = record.id.split('_')[0]
#         seq_target = str(record.seq).upper()
#         target_sequences.append((id_target, seq_target))
#         print(target_sequences)

# # For each target sequence, find the corresponding sequence in all.fasta
# # and extract the sequence to the right of the target sequence
# results = []
# for id_target, seq_target in target_sequences:
#     if id_target in all_sequences:
#         seq_all = all_sequences[id_target]
#         # Find where seq_target occurs in seq_all
#         index = seq_all.find(seq_target)
#         print(index, "Hiii")
#         if index != -1:
#             # Extract sequence to the right of the target sequence
#             start = index + len(seq_target)
#             right_seq = seq_all[start:]
#             results.append((id_target, right_seq))
#         else:
#             print(f"Target sequence not found in all.fasta sequence for ID {id_target}")
#             # Decide how to handle this case (empty sequence or skip)
#             results.append((id_target, ''))
#     else:
#         print(f"ID {id_target} not found in all.fasta")
#         # Decide how to handle this case (empty sequence or skip)
#         results.append((id_target, ''))

# # Write the results to a new FASTA file
# with open('output.fasta', 'w') as output_file:
#     for id_result, seq_result in results:
#         output_file.write(f">{id_target}\n")
#         output_file.write(f"{seq_result}\n")


########

# from Bio import SeqIO

# # Load all.fasta into a dictionary mapping IDs to records (including headers)
# all_sequences = {}
# with open('../dataset/test_structured_disordered.fasta', 'r') as all_file:
#     for record in SeqIO.parse(all_file, 'fasta'):
#         # Extract ID up to the first '|' character
#         id_all = record.id.split('|')[0]
#         all_sequences[id_all] = record  # Store the entire record

# # Load target.fasta into a list of tuples (ID, sequence)
# target_sequences = []
# with open('../dataset/test_L2_after_i2_cdhit.fasta', 'r') as target_file:
#     for record in SeqIO.parse(target_file, 'fasta'):
#         # Extract ID up to the first '_' character
#         id_target = record.id.split('_')[0]
#         # Convert the sequence to uppercase
#         seq_target = str(record.seq).upper()
#         target_sequences.append((id_target, seq_target))

# # For each target sequence, find the corresponding sequence in all.fasta
# # and extract the sequence to the right of the target sequence
# results = []
# for id_target, seq_target in target_sequences:
#     if id_target in all_sequences:
#         record_all = all_sequences[id_target]
#         seq_all = str(record_all.seq).upper()
#         # Find where seq_target occurs in seq_all
#         index = seq_all.find(seq_target)
#         if index != -1:
#             # Extract sequence to the right of the target sequence
#             start = index + len(seq_target)
#             right_seq = seq_all[start:]
#             # Keep the original header information
#             header = record_all.description
#             results.append((header, right_seq))
#         else:
#             print(f"Target sequence not found in all.fasta sequence for ID {id_target}")
#             # Decide how to handle this case (empty sequence or skip)
#             header = record_all.description
#             results.append((header, ''))
#     else:
#         print(f"ID {id_target} not found in all.fasta")
#         # Decide how to handle this case (empty sequence or skip)
#         header = id_target  # Use the ID as the header
#         results.append((header, ''))

# # Write the results to a new FASTA file, keeping the original headers
# with open('output.fasta', 'w') as output_file:
#     for header, seq_result in results:
#         output_file.write(f">{header}\n")  # Keep the original header info
#         output_file.write(f"{seq_result}\n")
