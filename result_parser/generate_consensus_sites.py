from Bio import AlignIO, SeqIO
from collections import Counter
import subprocess
import os
import argparse
import pandas as pd
import file_handlers as fh

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')
arguments.add_argument('-hr', '--hyphy_results',     help = 'Path to folder of hyphy results (including initial gene alignments)', required = True, type = str)

settings = arguments.parse_args()

current_directory = os.getcwd()
results_path = os.path.join(current_directory, settings.hyphy_results)
_PRE_MSA = "/data/shares/veg/dengue/hyphy-analyses/codon-msa/pre-msa.bf"
_POST_MSA = "/data/shares/veg/dengue/hyphy-analyses/codon-msa/post-msa.bf"


site_mappings_dir = os.path.join(results_path, "site_mappings")
unaligned_dir = os.path.join(site_mappings_dir, "unaligned")
pre_msa_dir = os.path.join(site_mappings_dir, "pre_msa")
mafft_dir = os.path.join(site_mappings_dir, "mafft")
post_msa_dir = os.path.join(site_mappings_dir, "post_msa")

for directory in [site_mappings_dir, unaligned_dir, pre_msa_dir, mafft_dir, post_msa_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)


def get_codon_consensus(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    consensus = []

    for i in range(0, alignment.get_alignment_length(), 3):  # Iterate codon-wise
        codons = [str(record.seq[i:i+3]) for record in alignment if '-' not in record.seq[i:i+3]]
        if codons:
            most_common_codon = Counter(codons).most_common(1)[0][0]
        else:
            most_common_codon = "---"  # Preserve gap if no consensus
        consensus.append(most_common_codon)

    return "".join(consensus)

def parse_fasta_biopython(file_path):
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}
    return sequences

def generate_CFEL_site_mapping(aligned_fasta):
    # Parse the aligned consensus sequences using Biopython
    aligned_consensus = parse_fasta_biopython(aligned_fasta)
    
    # Extract sequence names
    seq_names = list(aligned_consensus.keys())
    if len(seq_names) != 2:
        raise ValueError(f"Expected exactly two aligned consensus sequences (Processing file: {aligned_fasta})")
    # I know at least one of these serotypes had one gene with neither aligned sequences nor results, and that will cause problems here later
    
    seq1_name, seq2_name = seq_names
    seq1, seq2 = aligned_consensus[seq1_name], aligned_consensus[seq2_name]
    
    if len(seq1) % 3 != 0 or len(seq2) % 3 != 0:
        raise ValueError(f"Aligned sequences must have a length that is a multiple of 3 for codon-aware mapping (Processing file: {aligned_fasta})")
    if len(seq1) != len(seq2):
        raise ValueError(f"Aligned sequences must have the same length (Processing file: {aligned_fasta})")
    
    # Initialize site numbering
    site_num_1, site_num_2 = 0, 0
    mapping_info = []
    
    # Iterate through the aligned sequences codon-wise
    for i in range(0, len(seq1), 3):
        codon_seq1 = seq1[i:i+3]
        codon_seq2 = seq2[i:i+3]
        
        gap_in_seq1 = "-" in codon_seq1
        gap_in_seq2 = "-" in codon_seq2
        
        mapping_info.append([i // 3, site_num_1 if not gap_in_seq1 else "-", site_num_2 if not gap_in_seq2 else "-"])

        if not gap_in_seq1:
            site_num_1 += 1
        if not gap_in_seq2:
            site_num_2 += 1
    
    # Convert to DataFrame
    df_mapping = pd.DataFrame(mapping_info, columns=["Consensus_Site", f"{seq1_name}_Site", f"{seq2_name}_Site"])
    return df_mapping

def generate_BUSTED_site_mapping(aligned_fasta):
    # Parse the aligned consensus sequences using Biopython
    aligned_consensus = parse_fasta_biopython(aligned_fasta)
    
    # Extract sequence names
    seq_names = list(aligned_consensus.keys())

    # print(f"Processing file: {aligned_fasta}")

    # NS3 has no sequence alignment for DENV2, so skip DENV2
    if any("nonstructural_protein_NS3" in name for name in seq_names):
        if len(seq_names) != 3:
            raise ValueError("Expected exactly three aligned consensus sequences for nonstructural_protein_NS3")

        seq1_name, seq2_name, seq3_name = seq_names
        seq1, seq2, seq3 = aligned_consensus[seq1_name], aligned_consensus[seq2_name], aligned_consensus[seq3_name]
        
        if len(seq1) % 3 != 0 or len(seq2) % 3 != 0 or len(seq3) % 3 != 0:
            raise ValueError(f"Aligned sequences must have a length that is a multiple of 3 for codon-aware mapping (Processing file: {aligned_fasta})")
        if len(seq1) != len(seq2) or len(seq1) != len(seq3):
            raise ValueError(f"Aligned sequences must have the same length (Processing file: {aligned_fasta})")

        # Initialize site numbering
        site_num_1, site_num_2, site_num_3 = 0, 0, 0
        mapping_info = []
        
        # Iterate through the aligned sequences codon-wise
        for i in range(0, len(seq1), 3):
            codon_seq1 = seq1[i:i+3]
            codon_seq2 = seq2[i:i+3]
            codon_seq3 = seq3[i:i+3]
            
            gap_in_seq1 = "-" in codon_seq1
            gap_in_seq2 = "-" in codon_seq2
            gap_in_seq3 = "-" in codon_seq3
            
            mapping_info.append([
                i // 3, 
                site_num_1 if not gap_in_seq1 else "-", 
                site_num_2 if not gap_in_seq2 else "-", 
                site_num_3 if not gap_in_seq3 else "-"
            ])

            if not gap_in_seq1:
                site_num_1 += 1
            if not gap_in_seq2:
                site_num_2 += 1
            if not gap_in_seq3:
                site_num_3 += 1
        
        # Convert to DataFrame
        df_mapping = pd.DataFrame(mapping_info, columns=[
            "Consensus_Site", 
            f"{seq1_name}_Site", 
            f"{seq2_name}_Site", 
            f"{seq3_name}_Site", 
        ])
        return df_mapping
    else:
        if len(seq_names) != 4:
            raise ValueError(f"Expected exactly four aligned consensus sequences (Processing file: {aligned_fasta})")
        
        seq1_name, seq2_name, seq3_name, seq4_name = seq_names
        seq1, seq2, seq3, seq4 = aligned_consensus[seq1_name], aligned_consensus[seq2_name], aligned_consensus[seq3_name], aligned_consensus[seq4_name]
        
        if len(seq1) % 3 != 0 or len(seq2) % 3 != 0 or len(seq3) % 3 != 0 or len(seq4) % 3 != 0:
            raise ValueError(f"Aligned sequences must have a length that is a multiple of 3 for codon-aware mapping (Processing file: {aligned_fasta})")
        if len(seq1) != len(seq2) or len(seq1) != len(seq3) or len(seq1) != len(seq4):
            raise ValueError(f"Aligned sequences must have the same length (Processing file: {aligned_fasta})")
        
        # Initialize site numbering
        site_num_1, site_num_2, site_num_3, site_num_4 = 0, 0, 0, 0
        mapping_info = []
        
        # Iterate through the aligned sequences codon-wise
        for i in range(0, len(seq1), 3):
            codon_seq1 = seq1[i:i+3]
            codon_seq2 = seq2[i:i+3]
            codon_seq3 = seq3[i:i+3]
            codon_seq4 = seq4[i:i+3]
            
            gap_in_seq1 = "-" in codon_seq1
            gap_in_seq2 = "-" in codon_seq2
            gap_in_seq3 = "-" in codon_seq3
            gap_in_seq4 = "-" in codon_seq4
            
            mapping_info.append([
                i // 3, 
                site_num_1 if not gap_in_seq1 else "-", 
                site_num_2 if not gap_in_seq2 else "-", 
                site_num_3 if not gap_in_seq3 else "-", 
                site_num_4 if not gap_in_seq4 else "-"
            ])

            if not gap_in_seq1:
                site_num_1 += 1
            if not gap_in_seq2:
                site_num_2 += 1
            if not gap_in_seq3:
                site_num_3 += 1
            if not gap_in_seq4:
                site_num_4 += 1
        
        # Convert to DataFrame
        df_mapping = pd.DataFrame(mapping_info, columns=[
            "Consensus_Site", 
            f"{seq1_name}_Site", 
            f"{seq2_name}_Site", 
            f"{seq3_name}_Site", 
            f"{seq4_name}_Site"
        ])
        return df_mapping


for gene in fh.get_genes(results_path):
    concat_aln_path = os.path.join(results_path, "concat", "alignments", f"{gene}-nodups.fasta")
    with open(concat_aln_path, 'r') as concat_file:
        concat_aln = get_codon_consensus(concat_file)

    ### generate consensus sites between multi-clade and per-clade (for matching CFEL results to all other site results)
    clade_seqs = {}
    for clade in fh.all_clades:
        # if os.path.exists(os.path.join(site_mappings_dir, f"{gene}_{clade}_consensus_sites.tsv")):
        #     # no need to re-generate consensus sites if they already exist
        #     continue
        if gene == "nonstructural_protein_NS3" and clade == "DENV2":
            # (NS3 has no sequence alignment for DENV2, so skip it)
            continue

        clade_aln_path = os.path.join(results_path, clade, "alignments", f"{gene}-nodups.fasta")
        with open(clade_aln_path, 'r') as clade_file:
            clade_aln = get_codon_consensus(clade_file)
        clade_seqs[clade] = clade_aln


        CFEL_consensus_fasta = os.path.join(unaligned_dir, f"consensus_{gene}_{clade}.fasta")
        CFEL_aligned_consensus_fasta = os.path.join(post_msa_dir, f"aligned_consensus_{gene}_{clade}.fasta")

        # Save consensus sequences to FASTA format
        with open(CFEL_consensus_fasta, "w") as f:
            f.write(f">Consensus_{gene}_concat\n{concat_aln}\n")
            f.write(f">Consensus_{gene}_{clade}\n{clade_aln}\n")

        pre_msa_prot = os.path.join(pre_msa_dir, f'{gene}_{clade}_protein.fasta')
        pre_msa_nuc = os.path.join(pre_msa_dir, f'{gene}_{clade}_nuc.fasta')
        mafft_prot = os.path.join(mafft_dir, f'{gene}_{clade}.prot')
        # post_msa_out = os.path.join(post_msa_dir, f'{gene}_{clade}.msa')

        # run pre_msa script in hyphy to convert consensus nucleotide alignment to protein alignment
        try:
            print(f"Preprocessing gene sequences for {gene}, clade {clade}")
            subprocess.run([
            "hyphy", _PRE_MSA, 
            "--input", CFEL_consensus_fasta, 
            "--skip-realignment", "Yes",
            "--protein", pre_msa_prot, 
            "--rna", pre_msa_nuc
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running pre_msa for {gene}_{clade}_consensus_sites: {e}")

        # use mafft to align protein sequences
        try:
            print(f"Running MAFFT on {gene}, clade {clade}")
            subprocess.run(["mafft", 
                            "--quiet", 
                            pre_msa_prot], 
                            stdout=open(mafft_prot, "w"))
        except subprocess.CalledProcessError as e:
             print(f"Error running mafft for {gene}_{clade}_consensus_sites: {e}")

        # back-translate the protein alignment to nucleotide alignment
        # adds '_1' to the end of each fasta identifier for some godforsaken reason I don't care enough to debug
        try:
            print(f"Postprocessing gene sequences for {gene}, clade {clade}")
            subprocess.run([
            "hyphy", _POST_MSA, 
            "--protein-msa", mafft_prot,
            "--nucleotide-sequences", pre_msa_nuc, 
            "--output", CFEL_aligned_consensus_fasta
            ], check=True)
        except subprocess.CalledProcessError as e:
             print(f"Error running post_msa for {gene}_{clade}_consensus_sites: {e}")


        df_mapping = generate_CFEL_site_mapping(CFEL_aligned_consensus_fasta)

        # Save the mapping to a file
        df_mapping.to_csv(os.path.join(site_mappings_dir, f"{gene}_{clade}_consensus_sites.tsv"), sep="\t", index=False)

    ### generate consensus sites for all individual clades (for BUSTED site composition comparison across clades)
    # if os.path.exists(os.path.join(site_mappings_dir, f"{gene}_DENVall_consensus_sites.tsv")):
    #     # no need to re-generate consensus sites if they already exist
    #     continue

    # BUSTED_consensus_ref = os.path.join(site_mappings_dir, f"consensus_{gene}_{ref_clade}.fasta") # arbitrarily using DENV1 consensus as reference
    BUSTED_consensus_fasta = os.path.join(unaligned_dir, f"consensus_{gene}_DENVall.fasta")
    BUSTED_aligned_consensus_fasta = os.path.join(post_msa_dir, f"aligned_consensus_{gene}_DENVall.fasta")

    # Save consensus sequences to FASTA format
    with open(BUSTED_consensus_fasta, "w") as f:
        for c in fh.all_clades:
            if gene == "nonstructural_protein_NS3" and c == "DENV2":
                # (NS3 has no sequence alignment for DENV2, so skip it)
                continue
            f.write(f">Consensus_{gene}_{c}\n{clade_seqs[c]}\n")


    BUSTED_pre_msa_prot = os.path.join(pre_msa_dir, f'{gene}_DENVall_protein.fasta')
    BUSTED_pre_msa_nuc = os.path.join(pre_msa_dir, f'{gene}_DENVall_nuc.fasta')
    BUSTED_mafft_prot = os.path.join(mafft_dir, f'{gene}_DENVall.prot')

    # run pre_msa script in hyphy to convert consensus nucleotide alignment to protein alignment
    try:
        print(f"Preprocessing gene sequences for {gene}, DENVall")
        subprocess.run([
        "hyphy", _PRE_MSA, 
        "--input", BUSTED_consensus_fasta, 
        "--skip-realignment", "Yes",
        "--protein", BUSTED_pre_msa_prot, 
        "--rna", BUSTED_pre_msa_nuc
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running pre_msa for BUSTED {gene}_DENVall_consensus_sites: {e}")

    # use mafft to align protein sequences
    try:
        print(f"Running MAFFT on {gene}, DENVall")
        subprocess.run(["mafft", 
                        "--quiet", 
                        BUSTED_pre_msa_prot], 
                        stdout=open(BUSTED_mafft_prot, "w"))
    except subprocess.CalledProcessError as e:
            print(f"Error running mafft for BUSTED {gene}_DENVall_consensus_sites: {e}")

    # back-translate the protein alignment to nucleotide alignment
    try:
        print(f"Postprocessing gene sequences for {gene}, DENVall")
        subprocess.run([
        "hyphy", _POST_MSA, 
        "--protein-msa", BUSTED_mafft_prot,
        "--nucleotide-sequences", BUSTED_pre_msa_nuc, 
        "--output", BUSTED_aligned_consensus_fasta
        ], check=True)
    except subprocess.CalledProcessError as e:
            print(f"Error running post_msa for BUSTED {gene}_DENVall_consensus_sites: {e}")




    df_mapping = generate_BUSTED_site_mapping(BUSTED_aligned_consensus_fasta)

    # Save the mapping to a file
    df_mapping.to_csv(os.path.join(site_mappings_dir, f"{gene}_DENVall_consensus_sites.tsv"), sep="\t", index=False)


