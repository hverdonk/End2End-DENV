import sys
import re
import argparse

arguments = argparse.ArgumentParser(description='Get list of species that start with DENV*')
arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str)
arguments.add_argument('-o', '--output',  help = 'Directories ', required = True, type = str)
arguments.add_argument('-l', '--label',  help = 'Directories ', required = True, type = str)

settings = arguments.parse_args()
label = settings.label

with open(settings.input, 'r') as in_f:
	fasta_lines = in_f.readlines()

species_list = []
for ind, line in enumerate(fasta_lines):
	if line.startswith('>') and label in line:
		to_append = line.replace('>', '').rstrip()
		species_list.append(to_append)

species_string = ', '.join(species_list)
		

with open(settings.output, 'w') as out_f:
	out_f.write(species_string)

