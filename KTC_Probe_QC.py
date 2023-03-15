#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:07:29 2019

@author: u0121096
"""
from Bio import SeqIO
import time
import os
import argparse

Description = '''
'''

parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('probes_path', help='Insert full path to the file with the initial probe panel.')
parser.add_argument('refs_wo_IUPAC_path', help='Insert the full path to the file with references.')
parser.add_argument('initial_probes', help='Insert the full path to output directory.')
parser.add_argument('out_path', help='Insert the full path to output directory.')
parser.add_argument('-s','--similarity', help='Insert a value between 0 and 1 of the similarity threshold between probe and reference', required=True)

args = parser.parse_args()
probes_path = args.probes_path
refs_wo_IUPAC_path = args.refs_wo_IUPAC_path
initial_probes_path = args.initial_probes
reference_filename = os.path.splitext(os.path.basename(refs_wo_IUPAC_path))[0]
results_directory = args.out_path
log_path = os.path.join(results_directory, reference_filename + '_log.txt')
seq_similarity = args.similarity
Probe_design_preliminary = os.path.join(os.path.dirname(results_directory), 'Probe_design_preliminary')

start = time.time()

# Store the references in 'references'
references = []
for ref in SeqIO.parse(refs_wo_IUPAC_path, "fasta"):
	references.append(ref)
ref_len = len(references)

# Store probes in 'probes'
probes = []
for probe in SeqIO.parse(probes_path, "fasta"):
	probes.append(probe)
probe_length = len(str(probes[-1].seq))

#
#initial_probes = {}
#for probe in SeqIO.parse(initial_probes_path, "fasta"):
#	con = probe.id.split('_')
#	initial_probes.append(str(probe.seq))

# Function that score similarity of two sequences
def similarity(seq1, seq2):
	n_scoring = 0
	for n in range(len(seq1)):
		if seq1[n] == seq2[n]:
			n_scoring += 1
	return n_scoring/(len(seq1))




# Go through each reference and score the similarity of each probe to each possible position
# Write the probe, its best alignment coordinates, and alignment score to a file
for r in references:
	all_scores = []
	probes_final = []
	hyph_seq = ('-'*probe_length + r.seq.upper() + '-'*probe_length)
	for p_no, p in enumerate(probes):
		all_scores.append([])
		for n in range(len(r.seq)+probe_length+1):
			window = hyph_seq[n:n+probe_length]
			all_scores[p_no].append(similarity(p, window))
		with open(os.path.join(Probe_design_preliminary, r.id + '_probes_final.txt'), 'a') as out:
#			out.write('\n' + str(p.seq) + ',' + str(all_scores[p_no].index(max(all_scores[p_no]))-probe_length) + ',' + str(all_scores[p_no].index(max(all_scores[p_no]))) + ',' +  str(max(all_scores[p_no])) + ',' + '_'.join(p.id.split('_')[2:]))
			out.write('\n' + str(all_scores[p_no].index(max(all_scores[p_no]))-probe_length) + ',' + str(all_scores[p_no].index(max(all_scores[p_no]))) + ',' +  str(max(all_scores[p_no])) + ',' + '_'.join(p.id.split('_')[2:]))
	ref_end_time = time.time()
