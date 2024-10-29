import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import pandas as pd

# Load and parse the FASTA file
fasta_file = "L3_structured.fasta" # replace with the path to your fasta file
clades = []
lengths = []

for record in SeqIO.parse(fasta_file, "fasta"):
    # Extract clade from the description
    clade = record.description.split("_")[2]  # Adjust if clade position is different
    length = len(record.seq)
    
    clades.append(clade)
    lengths.append(length)

# Create DataFrame
data = pd.DataFrame({"Clade": clades, "Length": lengths})
data = data[data["Clade"] != "NA"]  # Exclude rows where clade is 'NA'

# Plot swarm plot with box plot
plt.figure(figsize=(12, 8))
sns.boxplot(data=data, x="Clade", y="Length", whis=[5, 95], fliersize=0, color="white")  # Box plot with neutral color
sns.swarmplot(data=data, x="Clade", y="Length", hue="Clade")  # Swarm plot overlay with hue
plt.xlabel("Clade", fontsize = 12)
plt.title(f"{fasta_file.split('.')[0]}")
plt.ylabel("Sequence Length", fontsize = 12)
# plt.title("Swarm and Box Plot of Sequence Length by Clade")
# plt.xticks(rotation=45)
plt.savefig(f"{fasta_file.split('.')[0]}.png")











# import matplotlib.pyplot as plt
# import seaborn as sns
# from Bio import SeqIO
# import pandas as pd

# # Load and parse the FASTA file
# fasta_file = "L2_after_i2_cdhit_disordered.fasta"  # replace with the path to your fasta file
# clades = []
# lengths = []

# for record in SeqIO.parse(fasta_file, "fasta"):
#     # Extract clade from the description
#     clade = record.description.split("_")[1]  # Adjust if clade position is different
#     length = len(record.seq)
    
#     clades.append(clade)
#     lengths.append(length)

# # Create DataFrame
# data = pd.DataFrame({"Clade": clades, "Length": lengths})

# # Plot swarm plot
# plt.figure(figsize=(10, 6))
# sns.swarmplot(data=data, x="Clade", y="Length")
# plt.xlabel("Clade")
# plt.ylabel("Sequence Length")
# plt.title("Swarm Plot of Sequence Length by Clade")
# plt.xticks(rotation=45)
# plt.savefig(f"{fasta_file}.png")
