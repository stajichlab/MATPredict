from Bio import Entrez, SeqIO

# Set your email address (required by NCBI)
Entrez.email = "kkell060@ucr.edu"

# Search query
query = '"Mating-type" OR "mating-type"'

# Maximum number of sequences to download
max_sequences = 50000

# Maximum sequence length
max_length = 1500

# Search for sequences
print("Searching for sequences...")
handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_sequences)
record = Entrez.read(handle)
handle.close()  # Close the handle after reading

# Get the list of sequence IDs
id_list = record["IdList"]

# Troubleshooting message
print(f"Found {len(id_list)} sequence IDs.")

myfile = open('MAT_loci.fasta', 'w')

# Set to keep track of species
species_set = set()

# Download the sequences
for seq_id in id_list:
    try:
        print(f"Fetching sequence with ID: {seq_id}")
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()  # Close the handle after reading
        
        # Extract species information
        species = record.annotations.get('organism', 'Unknown')
        
        # Check if the species is already in the set
        if species in species_set:
            print(f"Skipping sequence with ID: {seq_id} from species: {species}")
            continue
        
        # Add species to the set
        species_set.add(species)
        
        # Check the sequence length
        sequence_length = len(record.seq)
        if sequence_length > max_length:
            print(f"Sequence with ID: {seq_id} exceeds maximum length and was skipped.")
            continue

        # Convert the record to FASTA format
        sequence = record.format("fasta")
        myfile.write(sequence)
        
        # If the sequence is valid, print the success message
        print(f"Downloaded sequence with ID: {seq_id} from species: {species}")
    
    except Exception as e:
        print(f"Error processing sequence with ID: {seq_id}. Error: {e}")

myfile.close()
print("Finished downloading sequences.")
