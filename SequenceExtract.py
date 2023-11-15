from Bio import SeqIO
from Bio.Seq import Seq
import bisect
import warnings
import os
import sys

def get_YN(prompt):
    while True:
        user_input = input(prompt).lower()
        if user_input in ("y", "yes"):
            return True
        elif user_input in ("n", "no"):
            return False
        else:
            print("Invalid input. Please input 'y' for yes or 'n' for no.")

# Prompt the user for arguments
entrez_ids_file = input("Enter the path to the Entrez IDs file: ")
ref_genome_file = input("Enter the path to the reference genome file: ")
gtf_file = input("Enter the path to the GTF file: ")
output_fasta_file = input("Enter the path to the output FASTA file: ")
to_TTS = get_YN("Extract sequences to each gene's TTS (y/n): ")
upstream = int(input("Enter the number of upstream bases: "))
downstream = int(input("Enter the number of downstream bases (if to_TTS = y, this argument is ignored): "))
log_file = input("Enter the path to the log file: ")

# Create the log file if it doesn't exist
if os.path.exists(output_fasta_file):
    print(f"WARNING: {output_fasta_file} already exists.")
    sys.exit()
else:
    with open(output_fasta_file, 'w'):
        pass

# Create the log file if it doesn't exist
if os.path.exists(log_file):
    print(f"WARNING: {log_file} already exists.")
    sys.exit()

# Define all arguments in the header of the log file
else:
    with open(log_file, 'w') as log_file_handle:
        log_file_handle.write(f"Entrez ID file = {entrez_ids_file}\n")
        log_file_handle.write(f"Reference Genome file = {ref_genome_file}\n")
        log_file_handle.write(f"GTF file = {gtf_file}\n")
        log_file_handle.write(f"Ouptut FASTA file = {output_fasta_file}\n")
        log_file_handle.write(f"Extract sequences to each gene's TTS = {to_TTS}\n")
        log_file_handle.write(f"Upstream bases = {upstream}\n")
        log_file_handle.write(f"Downstream bases = {downstream}\n")
        log_file_handle.write(f"Log file = {log_file}\n\n")

# Load the reference genome
print("Loading reference genome...\n")
ref_genome = SeqIO.to_dict(SeqIO.parse(ref_genome_file, "fasta"))

# Read gene coordinates from GTF file
print("Parsing GTF file...\n")
gene_coordinates = {}
gene_end_lists = {}
with open(gtf_file, 'r') as gtf_file:
    for line in gtf_file:
        if line.startswith('#'):
            continue

        # Split the line by tab to get fields
        fields = line.strip().split('\t')

        # Check if feature is 'gene'
        feature_type = fields[2]
        if feature_type == 'gene':

            # Find the "db_xref" field
            db_xref_field = [field for field in fields[8].split(';') if "GeneID" in field]

            # Extract the gene ID
            if db_xref_field:
                gene_id = db_xref_field[0].split(':')[1].strip('"')
                gene_name = fields[8].split('"')[1]
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                gene_coordinates[gene_id] = (gene_id, gene_name, chrom, start, end, strand)

                # Add to gene end list
                if chrom not in gene_end_lists:
                    gene_end_lists[chrom] = [(end, gene_name)]
                else:
                    gene_end_lists[chrom].append((end, gene_name))

# Sort each gene end list
print("Sorting gene end lists...\n")
for chrom, end_list in gene_end_lists.items():
    gene_end_lists[chrom] = sorted(end_list, key=lambda x: x[0])

# Read Entrez Gene IDs from file
print("Checking Entrez IDs...\n")
with open(entrez_ids_file, 'r') as ids_file:
    entrez_ids = [line.strip() for line in ids_file]

# Extract sequences and write to output FASTA file
matched_genes = 0
trimmed_genes = 0
genes_not_found = []
print("Writing sequences to fasta file...\n")
with open(output_fasta_file, 'w') as output_fasta_file_handle, open(log_file, 'a') as log_file_handle:
    for gene_id in entrez_ids:
        if gene_id in gene_coordinates:
            matched_genes += 1
            gene_id, gene_name, chrom, start, end, strand = gene_coordinates[gene_id]

            # Calculate the upstream & downstream coordinates
            upstream_start = max(0, start - upstream)

            if to_TTS:
                downstream_stop = end
            else:
                downstream_stop = start + downstream

            # Check to see if upstream start collides with another gene
            upstream_index = bisect.bisect_left(gene_end_lists[chrom], (upstream_start,))
            start_index = bisect.bisect_left(gene_end_lists[chrom], (start,))
            if upstream_index != start_index:

                # Identify the collision gene
                collision_gene = gene_end_lists[chrom][start_index - 1]

                # Adjust the upstream start to the end of the collision gene
                upstream_start = collision_gene[0]

                # Print trimming report to the log file
                print(f"{gene_name} trimmed due to collision with {collision_gene[1]} at {collision_gene[0]}", file=log_file_handle)
                trimmed_genes += 1

            seq = ref_genome[chrom][upstream_start:downstream_stop]

            # Adjust sequence for reverse strand
            if strand == '-':
                seq = seq.reverse_complement()

            # Write sequence to fasta file
            seq.id = f"{gene_id}"
            seq.name = f"{gene_name}"
            seq.description = f"[organism=Homo sapiens][chromosome={chrom}] Gene={gene_name}, TSS={start}, coordinates={upstream_start}-{downstream_stop}, size={downstream_stop - upstream_start}"
            output_fasta_file_handle.write(seq.format("fasta"))

            # Write information to log file
            print(f"{seq}\n", file=log_file_handle)

        else:
            # Log information about genes not found in the GTF file
            genes_not_found.append(gene_id)
            print(f"Gene not found in GTF: {gene_id}")

    # Print gene IDs that were not found in the GFT file to the log file
    print("\n", file=log_file_handle)
    for gene_id in genes_not_found:
        print(f"ID not found in GTF: {gene_id}\n", file=log_file_handle)

    # Print the total number of trimmed genes to the log file
    print(f"\nTrimmed genes = {trimmed_genes}", file=log_file_handle)

# Print completion message
print(f"\nProcessing complete.\n{matched_genes} genes matched and written to:\n{output_fasta_file}")
