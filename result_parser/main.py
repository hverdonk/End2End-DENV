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
import os
from scipy.stats import chi2
import file_handlers as fh
import process_gene as pg
import traceback
import concurrent.futures


# Main script
if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Summarize hyphy analysis results for many genes into two csv files')
    arguments.add_argument('-hr', '--hyphy_results', help='Path to folder of hyphy results', required=True, type=str)
    settings = arguments.parse_args()

    current_directory = os.getcwd()
    results_path = os.path.join(current_directory, settings.hyphy_results)
    site_mappings_dir = os.path.join(results_path, "site_mappings")

    summary_fieldnames = ['gene', 'clade', 'N', 'T', 'dN/dS', 'sites', 'nt_conserved', 'aa_conserved',
                          'positive_sites', 'negative_sites', 'diff_sites',
                          'BUSTED_pval', 'BUSTED_omega3', 'BUSTED_prop_sites_in_omega3', 'RELAX_clade_K', 'RELAX_overall_pval']
    sites_fieldnames = ['gene', 'clade', 'site', 'consensus_site', 'composition', 'substitutions', 'majority_residue', 'diff_majority_residue',
                        'unique_aa', 'intensified_positive_selection', 'meme_marker', 'cfel_marker', 'prime_marker']

    genes = fh.get_genes(results_path)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(pg.process_gene, gene, results_path, current_directory, site_mappings_dir, summary_fieldnames, sites_fieldnames) for gene in genes]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                tb = traceback.format_exc()
                print(f'Generated an exception: {exc}\nTraceback: {tb}')


