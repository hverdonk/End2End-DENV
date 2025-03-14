"""
Combine analysis results

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)

Version:
    v0.0.1 (2021-01-17)
    v0.0.2 (2025-01-14)
    v0.0.3 (2021-02-26)


"""

import argparse
import csv
import random
import os
import json
import sys
import re
import math
import numpy
import gzip
import csv
from scipy.stats import chi2
from scipy.stats import mstats
import tree_helpers
import file_handlers as fh
from collections import Counter
import pandas as pd
import progressbar

bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)

random.seed ()

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-hr', '--hyphy_results',     help = 'Path to folder of hyphy results', required = True, type = str)
# arguments.add_argument('-o', '--output',     help = 'Output file prefix', required = True, type = str)
# arguments.add_argument('-fc', '--focal_clade',     help = 'Focal clade of interest, eg DENV1 (must exactly match label on tree)', required = True, type = str)
#arguments.add_argument('-oc', '--other_clade',     help = 'Clade to compare against the focal clade (e.g., background branches, must exactly match label on tree).', required = True, type = str)

settings = arguments.parse_args()

by_file = {}

timer = 0
count = 0
# tags = {}

current_directory = os.getcwd()
results_path = os.path.join(current_directory, settings.hyphy_results)
site_mappings_dir = os.path.join(results_path, "site_mappings")


def get_omega3 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 1)]




# i = 0
pv = 0.05

for gene in fh.get_genes(results_path):
  ####### Set up output files #######
  outfile_summary = os.path.join(current_directory, f"{gene}_summary.csv")
  outfile_sites = os.path.join(current_directory, f"{gene}_sites.csv")

  # Check if the files exist to determine whether to write headers
  # assert not os.path.exists(outfile_summary)
  # assert not os.path.exists(outfile_sites)

  summary_fieldnames = ['gene', 'clade', 'N', 'T', 'dN/dS', 'sites', 'nt_conserved', 'aa_conserved', 
                      'positive_sites', 'negative_sites', 'diff_sites', 
                      'BUSTED_pval', 'BUSTED_omega3', 'BUSTED_prop_sites_in_omega3', 'RELAX_clade_K', 'RELAX_overall_pval']
  sites_fieldnames = ['gene', 'clade', 'site', 'composition', 'substitutions', 'majority_residue', 'diff_majority_residue', 
                      'unique_aa', 'intensified_positive_selection','meme_marker', 'cfel_marker', 'prime_marker']


  with open(outfile_summary, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=summary_fieldnames)
    writer.writeheader()

  with open(outfile_sites, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=sites_fieldnames)
    writer.writeheader()
  ####### End set up output files #######


  # Load RELAX JSON
  try:
    relax_path = os.path.join(results_path, "concat", "RELAX", f"{gene}.RELAX.json")
    with open(relax_path, 'r') as relax_file:
      relax = json.load(relax_file)
  except FileNotFoundError:
    print(f"RELAX analysis not found for {gene}")
    relax = None

  # Load CFEL JSON
  cfel_path = os.path.join(results_path, "concat", "contrastFEL", f"{gene}.CFEL.json")
  with open(cfel_path, 'r') as cfel_file:
    cfel = json.load(cfel_file)


  ####### Preprocess CFEL data #######
  by_type = {}

  # get the tag for each branch (e.g., DENV1/DENV2/DENV3/DENV4 for the multi-serotype comparison)
  tested = cfel['tested']['0']  
  
  # group branch names by tag
  for b, t in tested.items():
      if not t in by_type:
          by_type[t] = []
      by_type[t].append (b)

  # Extract rows and headers from the cfel data
  data_rows = cfel["MLE"]["content"]["0"]
  headers   = cfel["MLE"]["headers"]

  # Build a quick lookup (short label -> index) for columns.
  header_map = {}
  for i, (short_label, _) in enumerate(headers):
      header_map[short_label] = i

  # Build lookups for each clade's beta and subs indices
  beta_idx_map = {}
  subs_idx_map = {}
  for c in fh.all_clades:
      beta_idx_map[c] = header_map[f"beta ({c})"]
      subs_idx_map[c] = header_map[f"subs ({c})"]
  ####### END Preprocess CFEL data #######


  for clade in fh.all_clades:
    gene_summary_dict = {'gene': gene, 'clade': clade}
    retag = re.compile (clade)
    invariant = [0,0]
    site_recorder = {}

    # Load results JSONs
    try:
      fel_path = os.path.join(results_path, clade, "FEL", f"{gene}.FEL.json")
      with open(fel_path, 'r') as fel_file:
        fel = json.load(fel_file)

      meme_path = os.path.join(results_path, clade, "MEME", f"{gene}.MEME.json")
      with open(meme_path, 'r') as meme_file:
        meme = json.load(meme_file)

      prime_path = os.path.join(results_path, clade, "PRIME", f"{gene}.PRIME.json")
      with open(prime_path, 'r') as prime_file:
        prime = json.load(prime_file)

      busted_path = os.path.join(results_path, clade, "BUSTED", f"{gene}.BUSTEDS.json")
      with open(busted_path, 'r') as busted_file:
        busted = json.load(busted_file)
    except FileNotFoundError:
      print(f"Analysis results not found for {gene} in clade {clade}")
      continue

    # Load consensus site mappings
    consensus_sites = os.path.join(site_mappings_dir, f"{gene}_{clade}_consensus_sites.tsv")
    consensus_sites_df = pd.read_csv(consensus_sites, sep='\t', dtype='Int64', na_values=['-', None])
    BUSTED_consensus_sites = os.path.join(site_mappings_dir, f"{gene}_DENVall_consensus_sites.tsv")
    BUSTED_consensus_sites_df = pd.read_csv(BUSTED_consensus_sites, sep='\t', dtype='Int64', na_values=['-', None])


    ####### begin RELAX analysis #######
    if relax:
      gene_summary_dict['RELAX_clade_K'] = relax["test results"]["relaxation or intensification parameter"][clade]
      gene_summary_dict['RELAX_overall_pval'] = relax["test results"]["p-value"]
    else:
      gene_summary_dict['RELAX_clade_K'] = "analysis not completed"
      gene_summary_dict['RELAX_overall_pval'] = "analysis not completed"
    ####### end RELAX analysis #######


    ####### begin CFEL analysis #######
    # For the current clade, print out summary for all branches associated with that clade
    # "Clade: N=number of branches for that clade, T=sum of branch lengths for the clade"
    print("*****")
    if clade in by_type:
        branches = by_type[clade]
        gene_summary_dict['N'] = len(branches)
        gene_summary_dict['T'] = sum([cfel['branch attributes']['0'][bn]['Global MG94xREV'] for bn in branches])
        print ("%s : N=%d, T=%g" % (clade, gene_summary_dict['N'], gene_summary_dict['T']))

        
    # Report dN/dS for the tested clade
    for k, r in cfel['fits']['Global MG94xREV']['Rate Distributions'].items():
      if clade in k:
        print ("%s: dN/dS = %g" % (k, r[0][0]))
        gene_summary_dict['dN/dS'] = r[0][0]


    # 1) Identify the index of key columns for the focal clade
    p_overall_idx     = header_map["P-value (overall)"]            # "P-value (overall)"
    beta_focal_label  = f"beta ({clade})"           # e.g. "beta (DENV1)"
    subs_focal_label  = f"subs ({clade})"           # e.g. "subs (DENV1)"
    beta_focal_idx    = header_map[beta_focal_label]
    subs_focal_idx    = header_map[subs_focal_label]

    # 2) Determine the other clades (hard-coded here; adjust as needed).
    other_clades = [c for c in fh.all_clades if c != clade]

    # 3) Build a lookup for pairwise p-value columns, e.g. "P-value for DENV1 vs DENV2"
    pairwise_pval_idx = {}
    for other in other_clades:
        # Two ways the column might appear:
        option_1 = f"P-value for {clade} vs {other}"
        option_2 = f"P-value for {other} vs {clade}"

        if option_1 in header_map:
            pairwise_pval_idx[other] = header_map[option_1]
        elif option_2 in header_map:
            pairwise_pval_idx[other] = header_map[option_2]
        # If neither is in the map, that pairwise p-value doesn't exist in the output

    # 4) Iterate over each site (row in data_rows).
    for k, row in enumerate(data_rows):
        # get consensus site number
        curr_site = fh.get_consensus_site_from_concat_site(consensus_sites_df, gene, k)

        # TODO: decide if I want to report sites that are invariant within a clade (using FEL) or across clades (as-is using CFEL)
        # mark if the site is invariant (ie first four entries in the list are 0)
        if max (row[0:4]) == 0.:
            invariant[0] += 1
        else:
            if max (row[1:4]) == 0.:
                invariant[1] += 1

        focal_beta = row[beta_focal_idx]
        overall_pval  = row[p_overall_idx]

        # -- PART 1: Check overall significance --
        # check if I need to know which beta is which / if I ever use the beta values again ----> (not in this script)
        # check if I need to know which subs is which / if I ever use the subs values again ----> (not from CFEL, not in this script)
        if overall_pval <= pv:
            # Prepare the list of [overall_pval, all betas, all subs]
            overall_info = [overall_pval]
            # Append beta values for each clade
            for c in fh.all_clades:
                overall_info.append(row[beta_idx_map[c]])
            # Append subs values for each clade
            for c in fh.all_clades:
                overall_info.append(row[subs_idx_map[c]])

            site_recorder.setdefault(curr_site, {})
            site_recorder[curr_site].setdefault("cfel", {})
            site_recorder[curr_site]["cfel"]["overall"] = {"site_info": overall_info}

            # If the focal clade's beta is larger than any other betas (test if selection was intensified in the focal clade compared to the other clades)
            if any(focal_beta > row[beta_idx_map[oc]] for oc in other_clades):
              site_recorder[curr_site]["cfel"]["overall"]["intensified"] = True
            else:
              site_recorder[curr_site]["cfel"]["overall"]["intensified"] = False

        # -- PART 2: Check pairwise significance (focal vs each other) --
        for other_clade in other_clades:
            if other_clade in pairwise_pval_idx:            
                pairwise_pval = row[pairwise_pval_idx[other_clade]]
                other_beta    = row[beta_idx_map[other_clade]]

                # Condition: pairwise p-value < pv and beta focal > beta other 
                if pairwise_pval <= pv: 
                    # e.g. "DENV1 vs DENV2": [pval, focal_beta, other_beta, focal_subs, other_subs]
                    comparison_key = f"{clade} vs {other_clade}"
                    comparison_info = [
                        pairwise_pval,
                        focal_beta,
                        other_beta,
                        row[subs_focal_idx],
                        row[subs_idx_map[other_clade]]
                    ]
                    
                    site_recorder.setdefault(curr_site, {})
                    site_recorder[curr_site].setdefault("cfel", {})
                    site_recorder[curr_site]["cfel"][comparison_key] = {"site_info": comparison_info} 

                    # is selection intensified in focal clade vs other clade?)
                    if focal_beta > other_beta:
                      site_recorder[curr_site]["cfel"][comparison_key]["intensified"] = True
                    else:
                      site_recorder[curr_site]["cfel"][comparison_key]["intensified"] = False


    gene_summary_dict['nt_conserved'] = invariant[0]
    gene_summary_dict['aa_conserved'] = invariant[1]
    print ("Sites %d (%d nt conserved, %d aa conserved)" % (cfel['input']['number of sites'],invariant[0],invariant[1]))
    ####### end CFEL analysis #######


    ####### begin FEL analysis #######
    gene_summary_dict['sites'] = fel['input']['number of sites']

    for i,k in enumerate(fel['MLE']['headers']):
        if k[0] == 'p-value': pvc = i
        if k[0] == 'alpha': ac = i
        if k[0] == 'beta': bc = i
            
    for k, i in enumerate(fel['MLE']['content']["0"]):
        # get consensus site number
        curr_site = fh.get_consensus_site_from_clade_site(consensus_sites_df, gene, clade, k)

        if i[pvc] <= pv:
            if not curr_site in site_recorder:
                site_recorder[curr_site] = {}
            site_recorder[curr_site]['fel'] = [i[pvc], i[bc] > i[ac]]
    ####### end FEL analysis #######
  

    ####### begin MEME analysis #######
    for i,k in enumerate(meme['MLE']['headers']):
        if k[0] == 'p-value': pvc = i
    
    for k, i in enumerate(meme['MLE']['content']["0"]):
        # get consensus site number
        curr_site = fh.get_consensus_site_from_clade_site(consensus_sites_df, gene, clade, k)

        if i[pvc] <= pv:
            if not curr_site in site_recorder:
                site_recorder[curr_site] = {}
            site_recorder[curr_site]['meme'] = [i[pvc]]
    ####### end MEME analysis #######


    ####### begin PRIME analysis #######
    prime_tags = []

    for i,k in enumerate(prime['MLE']['headers']):
        if k[0] == 'p-value': 
            pvc = i
            prime_tags.append ([0,'overall'])
        if k[0] == 'p1': 
            p1 = i
            prime_tags.append ([1, k[1].split (" ")[-1]])
        if k[0] == 'p2': 
            p2= i
            prime_tags.append ([2, k[1].split (" ")[-1]])
        if k[0] == 'p3': 
            p3 = i
            prime_tags.append ([3, k[1].split (" ")[-1]])
        if k[0] == 'p4': 
            p4 = i
            prime_tags.append ([4, k[1].split (" ")[-1]])
        if k[0] == 'p5': 
            p5 = i
            prime_tags.append ([5, k[1].split (" ")[-1]])
    
    prime_tags = [k[1] for k in sorted (prime_tags, key = lambda d: d[0])]
    
    for k, i in enumerate(prime['MLE']['content']["0"]):
        # get consensus site number
        # curr_site = consensus_sites_df.loc[consensus_sites_df[f"Consensus_{gene}_{clade}_Site"] == k, "Consensus_Site"].values[0]
        curr_site = fh.get_consensus_site_from_clade_site(consensus_sites_df, gene, clade, k)

        pvals = [i[j] for j in [pvc,p1,p2,p3,p4,p5]]
        if min (pvals) <= pv:
            if not curr_site in site_recorder:
                site_recorder[curr_site] = {}
            site_recorder[curr_site]['prime'] = [prime_tags[i] for i,k in enumerate (pvals) if k<= pv]
    ####### end PRIME analysis #######

    ####### begin BUSTED analysis #######
    pv = busted["test results"]["p-value"]
    if pv >= 0.05:
        print ("BUSTED p = %.2g" % (pv))
    else:
        print ("BUSTED p = %.2g, omega = %.2g (%.2g%%)" % (pv, get_omega3 (busted)['omega'], get_omega3 (busted)['proportion'] * 100.))

    gene_summary_dict['BUSTED_pval'] = pv
    gene_summary_dict['BUSTED_omega3'] = get_omega3 (busted)['omega']
    gene_summary_dict['BUSTED_prop_sites_in_omega3'] = get_omega3 (busted)['proportion'] * 100.
    ####### end BUSTED analysis #######


    print ("+    sites = %d" % (len([k for k in site_recorder.values() if 'meme' in k])))
    print ("-    sites = %d" % (len([k for k in site_recorder.values() if 'fel' in k and k['fel'][1] == False])))
    print ("diff sites = %d" % (len([k for k in site_recorder.values() if 'cfel' in k])))
    gene_summary_dict['positive_sites'] = len([k for k in site_recorder.values() if 'meme' in k])
    gene_summary_dict['negative_sites'] = len([k for k in site_recorder.values() if 'fel' in k and k['fel'][1] == False])
    gene_summary_dict['diff_sites'] = len([k for k in site_recorder.values() if 'cfel' in k])
    print("*****")

    # write gene summary statistics to output file
    with open(outfile_summary, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=summary_fieldnames)
        writer.writerow(gene_summary_dict)

    ####### begin SITES analysis #######
    def pr_composition (c):
        return " ".join (["%s/%d" % k for k in [k for k in sorted (list (c.items()), key = lambda d: -d[1]) if k[0] != '?']])

    def pr_subs (c):
        return " ".join (["%s/%d" % k for k in [k for k in sorted (list (c.items()), key = lambda d: -d[1]) if k[0] != '?']])
    
    # sort "sites" by site index, lowest first
    sites = sorted (([[i,k] for i,k in site_recorder.items() if 'meme' in k or 'cfel' in k]), key = lambda d: d[0])
    
    busted_internal_branches = {} # structure: {internal node: "test", leaf node: "background", ...}
    tree = tree_helpers.newick_parser (busted["input"]["trees"]["0"],{},busted_internal_branches)["json"]

    # print("\t\tmeme_marker\tcfel_marker\tprime_marker\tcomposition\tsubstitutions")

    for s in sites:
      other_clades = [c for c in fh.all_clades if c != clade]
      composition   = {}
      substitutions = {} 
      intensified_positive_selection = False
      unique_aa = None
      diff_majority_residue = False


      consensus_site_num = s[0]
      BUSTED_site_num = fh.get_clade_site_from_consensus_site(consensus_sites_df, gene, clade, consensus_site_num)
      # use BUSTED site num and BUSTED site index to get site nums for all other individual clades and feed them into the algorithm

      tree_helpers.traverse_tree (tree,  None, busted["substitutions"]["0"]["%d" % BUSTED_site_num], busted_internal_branches, composition, substitutions, clade)

      # get composition of other clades
      # consensus sites may not match between clades, since we're only doing pairwise alignments between each 
      # clade's consensus sequence and the concatenated clade consensus sequence
      # using consensus_site_num might screw things up for us here as a result
      for other_clade in other_clades:
        # other_consensus_sites = os.path.join(site_mappings_dir, f"{gene}_{other_clade}_consensus_sites.tsv")
        # other_consensus_sites_df = pd.read_csv(other_consensus_sites, sep='\t', dtype='Int64', na_values=['-', None])
        other_BUSTED_site_num = fh.get_clade1_site_from_clade2_site(BUSTED_consensus_sites_df, gene, clade, other_clade, BUSTED_site_num)
        #other_BUSTED_site_num = fh.get_clade_site_from_consensus_site(other_consensus_sites_df, gene, other_clade, consensus_site_num)

        other_busted_path = os.path.join(results_path, other_clade, "BUSTED", f"{gene}.BUSTEDS.json")
        with open(other_busted_path, 'r') as f:
            other_busted = json.load(f)

        other_internal_branches = {}
        other_tree = tree_helpers.newick_parser (other_busted["input"]["trees"]["0"],{},other_internal_branches)["json"]
        other_composition = {}

        if pd.isna(other_BUSTED_site_num):
            composition[other_clade] = {} # if the site is not present in the other clade, we'll just leave it empty
        else:
            tree_helpers.traverse_tree (other_tree, None, other_busted["substitutions"]["0"]["%d" % other_BUSTED_site_num], other_internal_branches, other_composition, {}, other_clade)
            composition[other_clade] = other_composition[other_clade]


      for k in fh.all_clades: 
        if k not in substitutions:
          substitutions[k] = {}   

      meme_marker = "%.3f" % s[1]['meme'][0] if 'meme' in s[1] else "-"
      
      if 'cfel' in s[1]:
        cfel_marker = ", ".join([f"{key}: {value['site_info'][0]:.3f}" for key, value in s[1]['cfel'].items()])
        intensified_selection = any([s[1]["cfel"][comparison_key]["intensified"] for comparison_key in s[1]["cfel"].keys()])
      else:
        cfel_marker = "-"
        intensified_selection = False

      prime_marker = ",".join (s[1]['prime']) if 'prime' in s[1] else "-"
      # a site is marked with * if the site is both positively selected and intensified in the focal clade relative to at least one other clade (beta_clade > beta_otherclade)
      if meme_marker != '-' and intensified_selection:
        site_marker = "*" 
        intensified_positive_selection = True
      else:
        site_marker = ""
      
      # a site is marked with § if there are amino-acid residues unique to the focal clade present
      # each BUSTEDS analysis is run on a single clade's alignment, NOT a multi-clade alignment
      # The genes can be different lengths in different clades, so "site 10" might be a totally different residue between clades
      # because the first residue got deleted in one clade and now residue 10 in clade 1 is residue 11 in clade 2, 
      # rather than there being different residues because of selection, which is what we're looking for
      # we'd have to run BUSTED on the multi-serotype alignments just to get the proper composition, and I'm not sure how practical it would be
      other_clade_compositions = set()
      for oc in other_clades:
          other_clade_compositions.update(composition[oc])
      if len(set(composition[clade]) - other_clade_compositions):
        site_marker += "§" 
        unique_aa = " ".join(set(composition[clade]) - other_clade_compositions)
      else:
        site_marker += ""

      # a site is marked with '#' if the majority residue is different between the focal clade and at least one other clade
      # the genes can be different lengths in different clades, so "site 10" might be a totally different residue just because they don't align, not because of selection
      focal_cons = sorted([[i1, i2] for i1, i2 in composition[clade].items()], key=lambda d: -d[1])[0]
      for oc in other_clades:
          other_cons = sorted([[i1, i2] for i1, i2 in composition[oc].items()], key=lambda d: -d[1])[0]
          if focal_cons[0] != other_cons[0]:
              site_marker += "#"
              diff_majority_residue = True
              break
      
    #   print ("Site\t%d %s\t%s\t%s\t%s\t%s\t%s" % (s[0]+1, 
    #                                               site_marker, 
    #                                               meme_marker , 
    #                                               cfel_marker, 
    #                                               prime_marker, 
    #                                               pr_composition (composition[clade]), 
    #                                               pr_subs (substitutions["test"]) if "test" in substitutions else "-"))
      
      gene_sites_dict = {
        'gene': gene,
        'clade': clade,
        'site': fh.get_clade_site_from_consensus_site(consensus_sites_df, gene, clade, consensus_site_num) + 1,
        'composition': pr_composition (composition[clade]),
        'substitutions': pr_subs(substitutions["test"]) if "test" in substitutions else "-",
        'majority_residue': focal_cons[0],
        'diff_majority_residue': diff_majority_residue,
        'unique_aa': unique_aa,
        'intensified_positive_selection': intensified_positive_selection,
        'meme_marker': meme_marker,
        'cfel_marker': cfel_marker,
        'prime_marker': prime_marker
      }

      with open(outfile_sites, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=sites_fieldnames)
        writer.writerow(gene_sites_dict)
    ####### end SITES analysis #######

