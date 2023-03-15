#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import os
from Bio import SeqIO
import argparse
import datetime

Description = '''
Synopsis
KTC_Consensus_generator.py is a script that generates a consensus based on an input of a file with the sequences to be aligned and condensed using muscle with '-gapopen -800'
It then calls the most frequent base (or hyphen) at each position.
A few notes:
\tThis script requires muscle to be installed and be able to be executed through the shell with a 'muscle' command
\tThis script was specifically designed to not consider or output IUPAC code. If a tie is found between two bases, a random base of the most frequent is called.
\tThis script prioritizes bases over hyphens. If, at a position, hyphens are the most numerous character, but the total count of all bases is higher, the most frequent base will be called. E.g:

\tA-A
\tA-A
\tA-A
\tAAA
\tAAA
\tATA
\tACA
\tAGA
\twill become AAA.

\tthis script has three optional arguments:
\t\t-nh: the consensus generated will not contain hyphens
\t\t-ac: the consensus generated is also saved together in an alignment with the input references.
\t\t-ka: the alignment generated from the references is saved

USAGE
This script is executed from the command line.
In the terminal, navigate to the folder with this script and then execute it. The format is as follows:
python KTC_Consensus_generator_muscle.py [path to one file with all sequences] [-nh, -ac, -ka]
	In order to insert the path to a file, drag and drop it in the terminal.

For additional help execute the command:
'python KTC_Sequence_Similarity.py -h'
The files generated are saved in the same folder as the input file.
\t'*_muscle.fasta' for the alignment including the consensus
\t'*_consensus.fasta' for the consensus sequence itself
'''


parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('references', help='Insert the full path to one file with (un)aligned references')
parser.add_argument('-nh', help='output a consensus that contains no hyphens', action="store_true")
parser.add_argument('-ac', help='saves an alignment for the references with the consensus', action="store_true")
parser.add_argument('-ka', help='keeps alignment of input references', action="store_true")
parser.add_argument('-l','--log_file', help='Insert the path to the log file')

args = parser.parse_args()
references_path = args.references
log_path = args.log_file
if args.nh:
	no_hyphens = True
else:
	no_hyphens = False
if args.ac:
	aligned_consensus = True
else:
	aligned_consensus = False
if args.ka:
	keep_alignment = True
else:
	keep_alignment = False

# These two functions are used to both print to the terminal and write to the log file. One of them (logger) includes the time of reporting.
def logger(msg):
	now = datetime.datetime.now()
	print(now.strftime("%Y-%m-%d %H:%M") + ' : ' + msg)
	with open(log_path, 'a') as out:
		out.write('\n' + now.strftime("%Y-%m-%d %H:%M") + ' : ' + msg)
def logger_no_clocker(msg):
	print(msg)
	with open(log_path, 'a') as out:
		out.write('\n' + msg)


#Defining the paths to the different files
alignment_path = os.path.splitext(references_path)[0] + '_muscle.fasta'
alignment_w_consensus = os.path.splitext(references_path)[0] + '_alignment_w_consensus.fasta'
consensus_path = os.path.splitext(references_path)[0] + '_consensus.fasta'

file_name = os.path.split(references_path)[1].split('.')[0]
results_directory = os.path.join(os.path.split(references_path)[0], 'Results') # Define the path to this new folder


# Sequences are loaded into a list as Biopython records
references = []
for ref in SeqIO.parse(references_path, 'fasta'):
	references.append(ref)

#The muscle command is run with twice the gapopen penalty
os.system("muscle -in %s -out %s -gapopen -800 -quiet" % (references_path,alignment_path))

#The aligned references are loaded into a list as Biopython records
aln_refs = []
for aligned_reference in SeqIO.parse(alignment_path, 'fasta'):
	aln_refs.append(aligned_reference)

# Non-degenerate IUPAC code is defined
bases = ['A', 'C', 'T', 'G', '-']

# Here we create a list with a number of lists equal to the length of the alignment
# Each item is, itself a list that will contain the count of each of the non-degenerate IUPAC code for that position
base_freqs = []
for n in range(len(aligned_reference.seq)):
	base_freqs.append([0, 0, 0, 0, 0])

# Degenerate IUPAC code is stored in a dictionary with those bases as keys for a list of its constituent bases
IUPAC = {
		'R': ['A', 'G'],
		'Y': ['C', 'T'],
		'S': ['G', 'C'],
		'W': ['A', 'T'],
		'K': ['G', 'T'],
		'M': ['A', 'C'],
		'B': ['C', 'G', 'T'],
		'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['A', 'C', 'G'],
		'N': ['A', 'G', 'C', 'T'],
		}

# =============================================================================
# Counting the bases
# =============================================================================
# We load the alignment and loop through each position for each reference.
# At each position we count the number of each of the non-degenerate bases
# Then we consider degenerate IUPAC code:
# If found, a score is added to each of the constituent bases.
# The score depends on the amount of characters that code covers. E.g. R = 0.5

for aln_ref in aln_refs:
	for n in range(len(aln_ref.seq)):
		for b in range(len(bases)):
			if aln_ref.seq[n].upper() == bases[b]:
				base_freqs[n][b] += 1
		for i in IUPAC:
			if aln_ref.seq[n].upper() == i:
				for b in range(len(bases)):
					for nt in IUPAC[i]:
						if bases[b] == nt:
							base_freqs[n][b] += 1/len(IUPAC[i])


# =============================================================================
# Generating the consensus
# =============================================================================
# The consenus is generated by looping through the score and assigning and appropriate base

consensus = "" # The consensus so far
ties = 0 # Variable that contains the number of ties between bases (that will be broken)
for x, n in enumerate(base_freqs):
	max_occurence = [i for i, x in enumerate(n) if x == max(n)] # max_occurence will be a list of items that are the most frequent
	if len(max_occurence) == 1 and sum(n[:-1]) >= n[-1]: # If there is only one most prevalent character and the sum of the count of bases is greater than the count of hyphens
		max_occurence = [i for i, x in enumerate(n[:-1]) if x == max(n[:-1])] # The most prevalent character is the most prevalent base (hyphens are excluded by excluding the last item of base_freqs)
		consensus += bases[random.choice(max_occurence)] # Add one of the most prevalent bases at random. In the case of a tie between hyphens and one base this is now just the one base
		if len(max_occurence) > 1: # But if the number of most frequent bases is still more than 1 (As would be the case for A, A, -, - , -, G, G)
			ties += 1
	elif len(max_occurence) == 1: # Otherwise, if there is only one most prevalent base (implying the sum of the frequencies of the bases is NOT greater than the count of hyphens): insert the most prevalent base (which is a hyphen)
		consensus += (bases[n.index(max(n))])
	elif 4 not in max_occurence: # Otherwise (if there are more than one most prevalent base) and a hyphen is NOT among them (index 4): choose one randomly from between the most prevalent bases
		consensus += bases[random.choice(max_occurence)]
		ties += 1
	else: # Otherwise (if there are more than one most prevalent base) and a hyphen is among them, choose one randomly from among the most prevalent bases
		consensus += bases[random.choice(max_occurence[:-1])]
		ties += 1

# Write a file with the alignment of the references together with the consensus if the use asked for it
if aligned_consensus == True:
	with open(alignment_path, 'r') as alignment:
		with open(alignment_w_consensus, 'a') as out:
			for line in alignment:
				out.write(line)
			out.write('>' + file_name + '_Consensus\n' + consensus)

# If the user did not ask for a consensus with hyphens write it without
# else write it without hyphens
if no_hyphens == False:
	with open(consensus_path, 'w') as out:
		out.write('>' + file_name + '_Consensus\n' + consensus)
else:
	with open(consensus_path, 'w') as out:
		out.write('>' + file_name + '_Consensus\n' + consensus.replace('-', ''))

# If the user asked to keep the alignment file
if keep_alignment == False:
	os.remove(alignment_path)

logger_no_clocker(str(len(references)) + ' references aligned to a length of: ' + str(len(aln_ref.seq)))
logger_no_clocker('Ties between most frequent nucleotides broken randomly: ' + str(ties))
if aligned_consensus == True:
	logger_no_clocker('Consensus in alignment written to: ' + alignment_w_consensus)
