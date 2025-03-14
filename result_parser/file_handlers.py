import os

all_clades = ["DENV1", "DENV2", "DENV3", "DENV4"]

def get_genes(path):
  genes = []
  # pick DENV1 FEL results arbitrarily, all clades should have the same genes
  for file in os.listdir(os.path.join(path, "DENV1", "FEL")):
    if file.endswith(".json"):
      gene = file.split(".")[0]
      genes.append(gene)
  return genes

# for matching CFEL site results to all other site results
def get_consensus_site_from_clade_site(consensus_sites_df, gene, clade, clade_site):
  try:
      return consensus_sites_df.loc[consensus_sites_df[f"Consensus_{gene}_{clade}_1_Site"] == clade_site, "Consensus_Site"].values[0]
  except IndexError:
      raise IndexError(f"index 0 is out of bounds for axis 0 with size 0 at {gene} {clade} {clade_site}")

def get_consensus_site_from_concat_site(consensus_sites_df, gene, clade, concat_site):
  try:
      return consensus_sites_df.loc[consensus_sites_df[f"Consensus_{gene}_concat_1_Site"] == concat_site, "Consensus_Site"].values[0]
  except IndexError:
      raise IndexError(f"index 0 is out of bounds for axis 0 with size 0 at {gene} concat (clade {clade}) {concat_site}")
  
def get_clade_site_from_consensus_site(consensus_sites_df, gene, clade, consensus_site_num):
  try:
      return consensus_sites_df.loc[consensus_sites_df["Consensus_Site"] == consensus_site_num, f"Consensus_{gene}_{clade}_1_Site"].values[0]
  except IndexError:
      raise IndexError(f"index 0 is out of bounds for axis 0 with size 0 at {gene} {clade} {consensus_site_num}")
  
# def get_concat_site_from_consensus_site(consensus_sites_df, gene, consensus_site_num):
#   try:
#       return consensus_sites_df.loc[consensus_sites_df["Consensus_Site"] == consensus_site_num, f"Consensus_{gene}_concat_Site"].values[0]
#   except IndexError:
#       raise IndexError(f"index 0 is out of bounds for axis 0 with size 0 at {gene} concat {consensus_site_num}")

# for matching BUSTED sites to do site composition comparisons
def get_clade1_site_from_clade2_site(consensus_sites_df, gene, clade1, clade2, clade1_site):
  try:
      return consensus_sites_df.loc[consensus_sites_df[f"Consensus_{gene}_{clade1}_1_Site"] == clade1_site, f"Consensus_{gene}_{clade2}_1_Site"].values[0]
  except IndexError:
      raise IndexError(f"index 0 is out of bounds for axis 0 with size 0 at {gene} {clade1} {clade2} {clade1_site} (BUSTED site matching)")


