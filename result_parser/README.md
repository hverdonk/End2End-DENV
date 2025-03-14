TO SUMMARIZE ALL HYPHY RESULTS FOR DENGUE, RUN AS FOLLOWS:

1. Build site indices between clades
`python3 generate_consensus_sites.py --hyphy_results clean_dataset_working_copy`

2. Summarize results for all hyphy analyses
`python3 result_summary.py --hyphy_results clean_dataset_working_copy`

Both programs expect the file structure output of snakefile-hyphy, 
snakefile-single-serotypes, snakefile-cfel-relax, and snakefile-pre to be unchanged 
(*i.e.*, one folder per serotype with separate subfolders for alignments and for each analysis result, 
and another concat folder with separate subfolders for alignments, for CFEL, and for RELAX results)
