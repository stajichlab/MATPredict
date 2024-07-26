from Bio import Entrez

# Set your email address (required by NCBI)
Entrez.email = "kkell060@ucr.edu"

# Search query
query = '"Mating-type" OR "mating-type"'

# Maximum number of sequences to download
max_sequences = 200

# Maximum sequence length
max_length = 1500

# Search for sequences
handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_sequences)
record = Entrez.read(handle)

# Get the list of sequence IDs
id_list = record["IdList"]

# Troubleshooting message
print(f"Found {len(id_list)} sequence IDs.")
myfile = open('MAT_loci.fasta', 'w')
# Download the sequences
for seq_id in id_list:
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    sequence = handle.read()
    print(f"Downloaded sequence with ID: {seq_id}")
    if len(sequence) < max_length:
        myfile.write(sequence)
    else:     
        print(f"Skipped sequence with ID: {seq_id}")
myfile.close()
# Close the handle
handle.close()