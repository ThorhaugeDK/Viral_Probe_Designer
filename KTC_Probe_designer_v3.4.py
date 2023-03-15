# -*- coding: utf-8 -*-

Description = '''
SYNOPSIS
KTC_Probe_designer is a script that generates a set of probes to reflect the diversity of a set of references. As input it needs two files:
  - A fasta file with references
  - A fasta file with a set of predefined, initial probes - these could be made, for instance, by splitting a consensus of your references up in smaller sequences of the desired probe length.
As probes may not contain IUPAC code, it is recommended that the initial probes be 'cleaned' of any IUPAC code within their sequences.
  Here I mean committing to a random base of one of the four bases contained within the degeneracy of the IUPAC code (i.e. ARA -> AAA / AGA).

Process
  1.The script 'cleans' all references of their IUPAC by converting a degenerate IUPAC code character into one of its constituent bases
  2.The script scans the first reference in the fasta file for the best position for each of the initial probes.
    - If a position for a probe is found that matches the reference sequence above or equal to a sequence similarity of a preset value (80%)
    - this part of the reference sequence is considered 'sufficiently covered' for the purposes of probe design and no probe is made for this position in the reference genome.
  3.The script then searches for gaps between 'sufficiently' covered regions
    - When a gap between these 'sufficiently covered' regions exceeds half the length of a probe:
      - This region will be bridged by a number of probes sufficient to cover the gap starting at the beginning of the gap
      - These probes are designed using the reference genome sequence and added to the initial probes.
    - The script then continues to the next reference using both the initial probes and the newly designed probes. In this fasion all references are looped over with an incresing amount of probes.
  4.All probes are matched again to each reference sequence for plotting their similarity and coverage.

  The output of the script is written to a folder in the directory of the reference file. This folder has the same name as the reference file but ends with '_probe_design'. This folder contains the following:
    - A fasta file with the reference sequences after they have been cleaned of their IUPAC code.
    - A fasta file with a list of all probes (initial as well as newly designed).
    - A png file for each reference with a graph for the similarity of each probe at its most similar position.

A few notes:
  - This script 'cleans' references in the same fashion prior to using them for probe design - All IUPAC code in references is converted into one of its constiuent bases at random
  - Be aware, that most probes are designed in the beginning. As more probes are designed, less are needed to cover the diversity of a subsequent reference.
    - It follows that references in the beginning of the fasta file will be used as a template for more probes than later references.
    - Consider placing references that you consider more general in their reflection of the diversity (such as consensus sequences) earlier than specific or rare sequences.

USAGE
This script is executed from the command line.
In the terminal, navigate to the folder with this script and then execute it. The format is as follows:
python KTC_Probe_designer.py [full path to one fasta file with all initial probes] [full path to one fasta file with all references] [optional: -s [sequence similarity threshold]
	In order to insert the path to a file, drag and drop it in the terminal.
'''

# =============================================================================
# Initial data management
# =============================================================================
# Script begins with importing modules need to capture elapsed time, manage filenames and directories and parse the command as executed in the command line.
import time
start = time.time()
import datetime
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('probes_path', help='Insert full path to the file with the initial probe panel.')
parser.add_argument('references_path', help='Insert the full path to the file with references.')
parser.add_argument('out_path', help='Insert the full path to output directory.')
parser.add_argument('-s','--similarity', help='Insert a value between 0 and 1 of the similarity threshold between probe and reference', required=True)
parser.add_argument('-t','--threads', help='Select the number of threads to be used', required=True)
parser.add_argument('-eg', '--exclude_genotype', nargs='*', help='list of space-seperated genotypes to be excluded for the purposes of probe design, but still tested against designed probes. If -cg is on, a consensus for these sequence this still be used to generate initial probes.')


args = parser.parse_args()
probes_path = args.probes_path
references_path = args.references_path
out_path = args.out_path
results_directory = out_path
seq_similarity = float(args.similarity)
threads = int(args.threads)
exclude_genotypes = args.exclude_genotype


reference_filename = os.path.splitext(os.path.basename(references_path))[0]
Probe_design_preliminary = os.path.split(references_path)[0]
path_to_plotter = './KTC_Probe_QC.py'

#The script checks to see whether the directories exists - if they do not, it creates them.
if not os.path.exists(results_directory):
	os.makedirs(results_directory)
if not os.path.exists(Probe_design_preliminary):
	os.makedirs(Probe_design_preliminary)

# Paths to the final probes new
final_probes_path = os.path.join(results_directory, reference_filename + '_final_probes.fasta')
refs_wo_IUPAC_path = os.path.join(results_directory, reference_filename + '_refs_wo_IUPAC.fasta')
log_path = os.path.join(results_directory, reference_filename + '_log.txt')

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

# The arguments used at launch are reported
logger_no_clocker('')
logger('\t-- Probe designer initialized --\n')
logger('\t-- Input --\n')
logger('Initial probes:\t ' + probes_path)
logger('References:\t\t ' + references_path)
logger('Similarity threshold: ' + str(seq_similarity))
logger('threads:\t\t ' + str(threads))
logger('Output directory:\t ' + results_directory)


#Probes and reference sequences are loaded into seperate lists, reference sequences are stored as BioPython 'SeqRecord' iterators. Initial probes are stored as SeqRecords. A copy of these is saved now as 'initial_probes'
references = []
probes = []
initial_probes = {}
initial_probe_count = 0
from Bio import SeqIO
for probe in SeqIO.parse(probes_path, "fasta"):
	probes.append(probe)
	con = '_'.join(probe.id.split('_')[2:])
	if con not in initial_probes:
		initial_probes[con] = []

	initial_probes[con].append(probe)
	initial_probe_count += 1

for ref in SeqIO.parse(references_path, "fasta"):
	references.append(ref)
#initial_probes = probes[:]


probe_length = len(str(probes[0].seq)) # Length of the probes are defined as the length of the first probe

# =============================================================================
# Cleaning IUPAC from reference sequences
# =============================================================================
# If launched with the probe_designer_master script, this step is redundant, but conducted anyway

logger_no_clocker('')
logger('\t-- Cleaning IUPAC from reference sequences (step 1/4) --\n')

#IUPAC code is defined in a dictionary with the keys containing lists of the constituent bases of the key.
refs_wo_IUPAC = {}
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

# Each reference is looped over with each character in its sequence that matches a character that is a key in the 'IUPAC' dictionary being replaced with a random base from the list of the dictionary with a key of that character.
ref_len = len(references)
ref_count = 0
refs_total_length = 0
import random
IUPAC_count = 0
for ref in references:
	ref_count += 1
	refs_wo_IUPAC[ref.id] = ''
	refs_total_length += len(ref.seq)
	for n in range(len(ref.seq)):
		if ref.seq.upper()[n] not in IUPAC:
			refs_wo_IUPAC[ref.id] += ref.seq.upper()[n]
		else:
			for I in IUPAC:
				if ref.seq.upper()[n] == I:
					refs_wo_IUPAC[ref.id] += random.choice(IUPAC[I])
					IUPAC_count += 1
	logger(ref.id + ' cleaned ' '(' + str(ref_count) + '/' + str(ref_len) + ' references)')

logger_no_clocker('')
logger('Cleaning of IUPAC completed - ' + str(IUPAC_count) + ' of ' + str(refs_total_length) + ' characters replaced (' + str(round((IUPAC_count/refs_total_length)*100, 2)) + '%)')

#The resulting 'cleaned' sequences are written to the output directory.
with open(refs_wo_IUPAC_path, 'w') as out:
		for ref in refs_wo_IUPAC:
			out.write('>' + ref + '\n' + refs_wo_IUPAC[ref] + '\n')

# The 'references' list that previously stored the raw reference sequences is emptied and then refilled with each of the 'cleaned' sequences reading from the newly created file.
references = []
for ref in SeqIO.parse(refs_wo_IUPAC_path, "fasta"):
	references.append(ref)


# =============================================================================
# Generating new probes
# =============================================================================

# This function is the core loop of probe design.
# It compares the sequence similarity of two input sequences.
# We are looping over this for each reference, for each probe for each position
# The scores are stored in the dictionary 'all_scores' under the key of the reference in question
def similarity(seq1, seq2):
	n_scoring = 0 # This score will increase every time the two sequences have identical bases in their position
	if len(seq1) != probe_length: # This is a safety measure. It should never happen that the length of the input sequence is not the length of the probe
		logger(probes[-1])
		logger('seq1: ' + seq1)
		logger('length of seq1: ' + str(len(seq1)))
		logger('Error: A probe was found not of the designated length. Terminating script.')
		quit()
	for n in range(len(seq1)):
		if seq1[n] == seq2[n]:
			n_scoring += 1
	return n_scoring/(len(seq1))

# This section stores the different variables necessary for the core loop of probe design
#all_scores = {} # This dictionary contains all scores for each probe, for each position, for each reference. Due to memory concerns, the information for one reference is removed before starting a new one
new_probes = {} # Contains the sequences for the new probes in a dictionary with references as the keys
refs_info = {} # Contains a list with the number of new probes and total probes made so far
probe_starts = [] # This list contains the base count at which a probe with >= seq_similarity starts covering the reference (remember all scores are relative to the extended reference (added hyphens at beginning and end))
ref_end_time = time.time() # The next three variables are for time keeping.
ref_start_time = time.time()
ref_duration_time = ref_end_time - ref_start_time

ref_len = len(references) # ref_len is the length of the variable 'references' ie: the number of references (If launched with the probe_designer_master this includes consensus sequences)

logger_no_clocker('')
logger('\t-- Starting design of new probes (step 2/4) --')
logger_no_clocker('')


# This loop is where probes are designed
# It identifies positions where probes align
# From this it decides where new probes are needed and makes them
for r_no, r in enumerate(references):
	if str(r.id)[0] in exclude_genotypes:
		logger(r.id + ' (' + str(r_no+1) + '/' + str(ref_len) + ' references) - ' + str(round((ref_end_time - start)/60, 2)) + ' minutes elapsed:')
		logger_no_clocker('\tReference skipped based on excluded genotypes\n')
		new_probes[r.id] = []
		refs_info[r.id] = [ref_duration_time, len(new_probes[r.id]), len(probes)]
	else:
		all_scores = [] # all_scores will contain a list with a number of lists equal to the number of possible positions a probe could align to the reference
		new_probes[r.id] = [] # A key with an empty list is generated for each reference 
		probe_log = [] # probe_log will contain the information on probes to be logger for the current reference
		hyph_seq = ('-'*probe_length + r.seq.upper() + '-'*probe_length) # hyphens are pre- and appended to the sequence. This allows the probe to 'slide in' when comparing similarities.
		for p_no, p in enumerate(probes):
			all_scores.append([])
		ref_start_time = time.time()
		probe_starts = [] # stores the starting position at which a probe aligns to a reference. '100' means from position 100 and forward (with a length of probe_length)
		probe_starts.append(0) # We put an ancor here saying the first probe_length characters in hyph_seq are covered (this is important for measuring the distances between covered regions)
		probe_starts.append(len(r.seq)+probe_length) # We put an anchor saying we are covering the last probe_length characters in hyph_seq (this is important for measuring the distances between covered regions)
		logger(r.id + ' (' + str(r_no+1) + '/' + str(ref_len) + ' references) - ' + str(round((ref_end_time - start)/60, 2)) + ' minutes elapsed:')
		for p_no, p in enumerate(probes): # For each probe we start a loop 'n' number of times. p_no will serve as the index of the list in the list all_scores that contain all scores for that probe at every position
			for n in range(len(r.seq)+probe_length+1): # n is the number of times we can slide along a length of the sequence length with a probe_lengths worth of hyphens at either end
				window = hyph_seq[n:n+probe_length] # window is a probe_length's worth of reference sequence starting from n'th base
				all_scores[p_no].append(similarity(str(p.seq).upper(), window)) # The similarity score for this reference for that probe, for that position is stored.
			if max(all_scores[p_no]) >= seq_similarity: # In the case that, that probe has a position for which its score is >= seq_similarity
				probe_starts.append(all_scores[p_no].index(max(all_scores[p_no]))) # Append that position to probe_starts (that position is the index of the highest scoring item relative to the hyphenated version of the reference)
		# With all the starts of well-aligning probes recorded we start identifying gaps in probe coverage for the current reference
		probe_starts.sort()
		for ps_no, ps in enumerate(probe_starts): # For each recorded start of a probe
			if ps != max(probe_starts): # Ignoring the highest number of probe starts (that is our anchor)
				probe_gap = min(x for x in probe_starts if x > ps) - ps # probe_gap is [the smallest number in the probe_starts variable that is larger than the current probe start] minus [the current probe start]. E.g. 360 - 240
				if probe_gap >= 1.5 * probe_length: # If probe_gap is >= 1.5*probe length (Note: This implies the stretch of uncovered sequence is >= 0.5*probe_length)
	
					probes_to_be_made = round(((probe_starts[ps_no+1] - probe_starts[ps_no])/probe_length)-1) # Probes_to_be_made determines how many whole probes are needed to bridge a gap.
					if probes_to_be_made >= 1:
						probe_log.append(str(probe_starts[ps_no]) + '\t' + str(probe_starts[ps_no+1]) + '\t' + str(probe_starts[ps_no+1] - probe_starts[ps_no] - probe_length) + '\t' + str(probes_to_be_made) + ' probe(s) have been designed')
						for n in range(probes_to_be_made): # For every probe that needs to be made
							if probe_starts[ps_no+1] == max(probe_starts) and n+1 == probes_to_be_made: # If we need to make a probe and the next probe_starts is the largest, make the probe from a probe_lengths worth of the last characters
								new_probes[r.id].append(str(r.seq.upper())[-probe_length:])
							else:
								new_probes[r.id].append(str(r.seq.upper())[ps+(probe_length*(n)):ps+(probe_length*(n+1))]) # Make a set of n probes end to end starting from the current probe start


		if not len(probe_log) == 0: # If new probes had to be made for this reference display the gaps in coverage
			logger_no_clocker('\nstart\tend\tgap')
			for gap in probe_log:
				logger_no_clocker(gap)
	
		if len(new_probes[r.id]) != 0:
			for p_no, p in enumerate(new_probes[r.id]): # Add the probes that were made for the reference to the list of all probes
				with open(os.path.join(Probe_design_preliminary, r.id + '_new_probes.fasta'), 'a') as new_probes_fasta:
					new_probes_fasta.write('\n>Probe_' + str(len(probes)+p_no+1) + '_' + r.id + '\n' + p)
	
			for probe in SeqIO.parse(os.path.join(Probe_design_preliminary, r.id + '_new_probes.fasta'), 'fasta'):
				probes.append(probe)
	
	#		probes.append(p)
		ref_end_time = time.time()
		ref_duration_time = ref_end_time -ref_start_time
		refs_info[r.id] = [ref_duration_time, len(new_probes[r.id]), len(probes)]
		# Log how many probes aligned to the reference and how many new probes had to be made
		logger_no_clocker('\n' + str(len(probe_starts)-2) + ' probes aligned to the reference.')
		logger_no_clocker(str(len(new_probes[r.id])) + ' new probes designed to accomodate this reference. New probe total: ' + str(len(probes)))
		logger_no_clocker('')


#When all probes have been made, write them to a new fasta file with their IDs

#print(probes)
probe_no = 0
with open(final_probes_path, 'w') as out:
	for p in probes:
		probe_no += 1
#		print(p)
		out.write('>Probe_' + str(probe_no) + '_' + '_'.join(p.id.split('_')[2:]) + '\n' + str(p.seq) + '\n')

# =============================================================================
# Analyzing the results
# =============================================================================
logger('\t-- Starting analysis of all probes (step 3/4) --\n')


# We split the references into a number of buckets depending on the number of references and threads to be engaged
import math
import subprocess

ref_buckets = {} # ref_buckets will contain reference IDs with the bucket number as their keys
bucket_total = math.ceil(ref_len/threads) # this will be the number of buckets. E.g. 10 references using 3 threads = 3.33 => 4 buckets

# Create n number of buckets (keys in ref_buckets)
for n in range(1, bucket_total+1):
	ref_buckets[n] = []

# References are split into buckets
bucket_no = 1 
for ref in references:
	SeqIO.write(ref, os.path.join(Probe_design_preliminary, ref.id + '.fasta'), 'fasta') # References are written individually to new files
	if bucket_no == len(ref_buckets): # If the current bucket_no is the number of buckets, write that reference to that bucket and reset bucket_no to 1
		ref_buckets[bucket_no].append(ref)
		bucket_no = 1
	else: # Else write the reference to the current bucket_no and increase the bucket_no by 1
		ref_buckets[bucket_no].append(ref)
		bucket_no += 1

logger('References split into ' + str(bucket_total) + ' buckets.')

# A command is written for each bucket in the format: command1 &\n command2 &\n... wait
# Doing this means engaging a number of threads equal to the commands written, waiting until all have completed for going to next bucket
for bucket_no, bucket in enumerate(ref_buckets):
	logger_no_clocker('')
	logger('Now aligning probes against references in bucket ' + str(bucket_no+1) + '/' + str(bucket_total))
	mp_command = ''
	for ref in ref_buckets[bucket]:
		logger_no_clocker('\t' + ref.id)
		mp_command += 'python ' + path_to_plotter + ' ' + final_probes_path + ' ' + os.path.join(Probe_design_preliminary, ref.id + '.fasta') + ' ' + probes_path + ' ' + results_directory + ' -s ' + str(seq_similarity) +  ' &\n'
	mp_command +='wait'
	subprocess.run(mp_command, shell=True)


logger_no_clocker('')

# =============================================================================
# Drawing plots of probe coverage and similarity
# =============================================================================
logger('\t-- Drawing plots of probe coverage (step 4/4) --')

# Priscillas colours
#colour_dict = {
#		'con_of_cons': 'black',
#		'1': '#9ca7ff',
#		'2': '#fae173',
#		'3': '#c9fff3',
#		'4': '#ffa99c',
#		'5': '#1da384',
#		'6': '#a14d82',
#		'7': 'magenta',
#		'8': '#ffffff'
#		}


colour_dict = {
		'consensus': 'black',
		'con_of_cons': 'white',
		'1': '#4f4dff', # blue
		'2': '#9d4dff', # purple
		'3': '#ff6666', # red
		'4': '#fff766', # yellow
		'5': '#a1ff66', # green
		'6': '#66ffed', # teal
		'7': '#ffc266', # orange
		'8': '#ff66cc' # pink
		}
	
# Create a graph for each probe, reading the necessary information from the _probes_final.txt files
for ref in references:
	logger_no_clocker('\t' + ref.id)
#	try:
	covered_positions = ['r']*len(ref.seq)
	with open(os.path.join(Probe_design_preliminary, ref.id + '_probes_final.txt'), 'r') as probes_final:
		fig, ax= plt.subplots()
		fig.set_size_inches(30,10)
		ax.set_title(ref.id, size = 40, y=1.08)
		genotypes_seen = []



		for p_no, line in enumerate(probes_final):

			if not line.startswith('\n'):
				x_coord = int(line.split(',')[0])
				y_coord = float(line.split(',')[2])
#				print(x_coord, y_coord)
				if y_coord >= seq_similarity:
#					print(y_coord)
					for pos in range(x_coord, x_coord+probe_length):
						if pos < len(covered_positions) and pos >= 0:
							covered_positions[pos] = 'g'

			if not line.startswith('\n'):
				x = [int(line.split(',')[0]), int(line.split(',')[1])]
				y = [float(line.split(',')[2]), float(line.split(',')[2])]
				probe_source = line.split(',')[3].replace('\n', '')
				if probe_source in initial_probes:
					plt.plot(x,y, marker="|", color=colour_dict[probe_source], label = probe_source if probe_source not in genotypes_seen else "")
					if probe_source not in genotypes_seen:
						genotypes_seen.append(probe_source)
				else:
					plt.plot(x,y, marker="|", color='grey', label = 'input reference' if 'input reference' not in genotypes_seen else "")
					if 'input_reference' not in genotypes_seen:
						genotypes_seen.append('input reference')


#		print(covered_positions)
		for index, colour in enumerate(covered_positions):
			if colour == 'g':
				plt.axvline(x = index, color = 'green', alpha = 0.02)
			else:
				plt.axvline(x = index, color = 'red', alpha = 0.03)




		plt.xticks(np.arange(0, len(ref.seq), 1000), size = 20)
		plt.yticks(np.arange(0, 1.05, 0.05), size = 10)
		plt.axvline(x = 0, color = 'black')
		plt.axvline(x = len(ref.seq), color = 'black')
		plt.axhline(y = seq_similarity, color = 'red', linestyle='dashed', linewidth=0.5)
		ax.set_xlabel("position in reference", size = 30)
		ax.set_ylabel("probe-reference similarity", size = 30)
		plt.legend(loc='lower right', )
		plt.tight_layout()
		fig.savefig(os.path.join(results_directory, ref.id.replace('/', '_') + '.png'))
		plt.close()
#	except:
#		logger_no_clocker(ref.id + ' failed.')


# Create a new graph - its size depends on the number of references
fig, ax = plt.subplots(figsize=(3+ref_len*0.2,8))

x = np.array([x for x in range(ref_len+1)]) # A list of numbers from 0
x_ticks = [] # This variable holds the strings of the reference names in a list
new_p = [] # This variable holds the number of new probes in a list
total_p = [] # This variable holds the accumulated number of probes in a list
for ref in refs_info: # refs_info holds information on the number of new probes per reference
	new_p.append(refs_info[ref][1])
	x_ticks.append(ref)
	total_p.append(refs_info[ref][2])
x_ticks = np.array(['Initial probes'] + x_ticks)

total_p = np.array([initial_probe_count] + total_p)
for ref_no, ref in enumerate(refs_info):
	plt.xticks(x, x_ticks, rotation=45, ha='right')
	plt.plot(x,total_p, color='#ff8c00')
	if ref_no == ref_len-1:
		plt.text(x[-1]+0.1, total_p[-1], total_p[-1], fontsize=20)

ax.set_ylabel('Total probes', color='#ff8c00', fontsize=20)
plt.ylim(0, probe_no+10)
ax.xaxis.grid(True, alpha=0.6)

ax2 = ax.twinx()
min_new_p = min(new_p)
max_new_p = max(new_p)

new_ps = np.array([np.NaN] + new_p)
for ref in refs_info:
	plt.plot(x,new_ps, color='steelblue')
plt.yticks(np.arange(0, max_new_p+1, 1))
plt.xticks(x, x_ticks, rotation=45, ha='right')
ax2.set_ylabel('New probes', color='steelblue', fontsize=20)
plt.grid(True, alpha=0.6)
ax.set_title(os.path.basename(results_directory), size = 40, y=1.08)
fig.tight_layout()
fig.savefig(os.path.join(results_directory, 'probe_count.png'), dpi=300)
plt.close()



from statistics import median

logger_no_clocker('')
logger('Now testing probe-to-probe similarity')
all_similarities = []
max_similarities = []

# Probe to probe similarity is checked using the similarity function.
for p_no, p in enumerate(probes):
	p_excised = probes[:p_no] + probes[p_no+1:] # p_excised is a list of all probes without the probe being checked for that loop
	all_similarities = []
	for pe_no, pe in enumerate(p_excised):
		for n in range(probe_length+probe_length+1):
			window = ('-'*probe_length + pe.upper() + '-'*probe_length)[n:n+probe_length]
			all_similarities.append(similarity(p.upper(), window))
	max_similarities.append(max(all_similarities))

import matplotlib.pyplot as plt
median = median(max_similarities)
max_si = max(max_similarities)
min_si = min(max_similarities)

# Creating the boxplot of probe-probe similarities
plt.figure(figsize = (5, 7))
plt.boxplot(max_similarities)
plt.yticks(np.arange(0, 1.1, 0.1))
plt.axhline(y = seq_similarity, color = 'red', linewidth=0.5, linestyle='dashed')
plt.title('Probe - Probe similarity', fontsize = 15)
plt.ylabel('Similarity', fontsize = 15)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.text(1.1, median, round(median, 2))
plt.text(1.1, max_si, round(max_si, 2))
plt.text(1.1, min_si, round(min_si, 2))
plt.savefig(os.path.join(results_directory, 'probe_QC_boxplot.png'), dpi=300)
plt.close()


end = time.time()
logger('Execution time: ' + str(round((end-start)/60, 2)) + ' minutes.')

