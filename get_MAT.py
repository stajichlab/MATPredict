from Bio import Entrez, SeqIO

# Set your email address (required by NCBI)
Entrez.email = "kkell060@ucr.edu"

# Search query
query = '"Mating-type" OR "mating-type"'

# Maximum number of sequences to download
max_sequences = 20000

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

# Set to keep track of species
species_set = set()

# Download the sequences
for seq_id in id_list:
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    
    # Extract species information
    species = record.annotations.get('organism', 'Unknown')
    
    # Check if the species is already in the set
    if species in species_set:
        print(f"Skipping sequence with ID: {seq_id} from species: {species}")
        continue
    
    # Add species to the set
    species_set.add(species)
    
    # Convert the record to FASTA format
    sequence = record.format("fasta")
    
    # Truncate the sequence if it exceeds the maximum length
    header, seq = sequence.split('\n', 1)
    seq = seq.replace('\n', '')
    if len(seq) > max_length:
        seq = seq[:max_length]
    
    # Write the header and truncated sequence to the file
    myfile.write(f"{header}\n")
    myfile.write('\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n')
    
    print(f"Downloaded sequence with ID: {seq_id} from species: {species}")

myfile.close()
# Close the handle
handle.close()