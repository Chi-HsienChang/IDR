import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import pandas as pd

# Define file paths
fasta_file_L3 = "L3_IDRs.fasta"  # Replace with the actual path if different
fasta_file_L2 = "L2_IDRs.fasta"  # Replace with the actual path if different

# Function to parse FASTA and extract clade and sequence length
def parse_fasta(file, label):
    clades = []
    lengths = []
    labels = []
    for record in SeqIO.parse(file, "fasta"):
        clade = record.description.split("_")[2]
        length = len(record.seq)
        clades.append(clade)
        lengths.append(length)
        labels.append(label)
    return pd.DataFrame({"Clade": clades, "Length": lengths, "File": labels})

# Parse each file and combine into a single DataFrame
data_L3 = parse_fasta(fasta_file_L3, "L3")
data_L2 = parse_fasta(fasta_file_L2, "L2")
data = pd.concat([data_L3, data_L2])

# Exclude rows where clade is 'NA'
data = data[data["Clade"] != "NA"]

# Plot swarm and box plot in the same figure
plt.figure(figsize=(14, 8))
sns.boxplot(data=data, x="Clade", y="Length", hue="File", whis=[5, 95], fliersize=0, color="white", linewidth=1.5)
sns.swarmplot(data=data, x="Clade", y="Length", hue="File", dodge=True, palette="Set2", marker="o")

# Customize plot
plt.xlabel("Clade")
plt.ylabel("Sequence Length")
plt.title("Swarm and Box Plot of Sequence Length by Clade for L3 and L2")
plt.xticks(rotation=45)
plt.legend(title="File", bbox_to_anchor=(1.05, 1), loc='upper left')

# Save the figure
plt.savefig("IDRs.png")
# plt.show()
