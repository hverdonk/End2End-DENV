from Bio import SeqIO
import sys
import argparse
import os

arguments = argparse.ArgumentParser(description = 'Get seq from genbank file')
arguments.add_argument('-i', '--input_file', help='Genbank file (full path)', required=True, type=str)
arguments.add_argument('-o', '--output_dir', help='Fasta file output (full path)', required=True, type=str)
arguments.add_argument('-c', '--csv_file', help='Txt file with protein names, start and stop', required=True, type=str)
#arguments.add_argument('-e', '--end', help='CDS end', required=True, type=str)
settings = arguments.parse_args()

def extract_fasta(genbank_file, start, end, product, output_fasta):
    """
    Extracts a specified sequence from a GenBank file and writes it to a FASTA file.
    
    Args:
        genbank_file (str): Path to the input GenBank file.
        start (int): Start coordinate (1-based inclusive).
        end (int): End coordinate (1-based inclusive).
        product (str): Product name to filter features (e.g., "capsid protein C").
        output_fasta (str): Path to the output FASTA file.
    """
    # Parse the GenBank file
    record = SeqIO.read(genbank_file, "genbank")
    
    # Check for features matching the product
    for feature in record.features:
#        if feature.type == "CDS" and "product" in feature.qualifiers:
#            if product in feature.qualifiers["product"]:
                # Extract the sequence
         extracted_seq = record.seq[start - 1:end]  # Convert to 0-based indexing
                
                # Write to FASTA
         with open(output_fasta, "w") as fasta_file:
             fasta_file.write(f">{record.id}|{product}|{start}-{end}\n")
             fasta_file.write(str(extracted_seq) + "\n")
         
         print(f"Sequence extracted and saved to {output_fasta}")
         return
    
    print(f"No matching feature with product '{product}' found in {genbank_file}.")

# Example usage
if __name__ == "__main__":
    # Provide input details here or via command line arguments
    genbank_file = settings.input_file #"dengue-virus-1.gb"
    with open(settings.csv_file, 'r') as csv_file:
        lines = csv_file.readlines()
    for line in lines:
        line = line.rstrip()
        line = line.split(',')
        output_fasta =  os.path.join(settings.output_dir, line[0] + '.fasta') #settings.output_file #"capsid_protein_C.fasta"
        start = int(line[1]) #int(settings.start) #95
        end = int(line[2]) #int(settings.end) #394
        product = line[0] #"capsid protein C"
        sanitized_product = re.sub(r'[\/:*?"<>|]', '_', product) # "capsid_protein_C"
        # Run the extraction
        extract_fasta(genbank_file, start, end, sanitized_product, output_fasta)
