#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import os
import subprocess
import textwrap
import datetime
import sys
import argparse

Description = '''
SYNOPISIS
'KTC_probe_designer_master_v1.1.py' is a fully automatic script that generates a set of probes that reflect the diversity of a set of reference sequences.
As input, KTC_probe_designer_master_v1.3.py requires only a set of references (which have headers designating their subtype (ie '>1a_EU123456' for a subtype 1a)) and a decimal number between 0 and 1 for probe to reference similarity that defines when a part of the reference is dissimilar enough from all probes that it justifies making a new probe.
The script then runs through the following steps in a sequential manner:
  - All ambigious IUPAC code is broken down to one of its constituent bases randomly to comply with probe design limitations
  - The 3' poly(T) tract is excised by scanning the last 400 nt of each sequence, cutting when two 'T's in a row AND 13 of the following 20 bases are 'T's
    - The excision point is reported for manual confirmation
  - Sequences are separated by their subtype designated by their headers and written into distinct fasta files
  - The script: 'KTC_Consensus_Generator_Muscle.py' generates one consensus for each subtype:
    - The fasta files each containing all references of a specific subtype are aligned using muscle with a -gapopen penalty of -800
    - A consensus is generated using the most frequent base
  - All subtype consensuses are written to the same fasta file
  - A consensus of all the consensuses is generated using the script 'KTC_Consensus_Generator_Muscle.py'
  - The resulting sequence is split into words of 120 nt. The last word may not contain 120 nt, so it is removed and replaced by the last 120 characters of the sequence. These ~79 sequences constitute our initial probes.
  - The script 'KTC_Probe_designer_v3.2d.py' is executed using the input references, the generated consensuses and the similarity threshold determined when running this script (see 'KTC_Probe_designer_v3.2d.py' for details)

The output is written to the same directory as the input references and written into two subdirectories:
  - 'Probe_designer_preliminary': for all preparation work prior to the actual probe generation
  - 'Results': for all figures, the log file, the references with IUPAC converted, and the designed probes themselves (Results/[reference file name]_new_probes.txt)

USAGE
This script is executed from the command line.
In the terminal, navigate to the folder with this script and then execute it. The format is as follows:
python KTC_probe_designer_master_v1.1.py [full path to one fasta file with all references] [-s [similarity threshold value as a decimal between 0 and 1 (use '.') as separator] (optional, default is 0.8)] [-t [number of threads to be engaged] (optional, default is 1)]
	In order to insert the path to a file, drag and drop it in the terminal.
'''



# =============================================================================
# Initialization and housekeeping
# =============================================================================
# argparse permits the script to be run from the terminal. The path to the references is defined as references and the input similarity threshold as the float of the input
parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('references_path', help='Insert the full path to the file with references.')
parser.add_argument('-s','--similarity', help='Insert a value between 0 and 1 of the similarity threshold between probe and reference', default=0.80)
parser.add_argument('-t','--threads', help='Select the number of threads to be used', default=1)
parser.add_argument('-c', '--consensus', help='The consensus sequence of all samples is used to generate initial probes', action='store_true')
parser.add_argument('-cc', '--consensus_of_consensus', help='The consensus of all subtype consensuses is used to generate initial probes', action='store_true')
parser.add_argument('-cg', '--consensus_of_genotypes', help='Each consensus of a genotype is used to generate initial probes', action='store_true')
parser.add_argument('-cs', '--consensus_of_subtypes', help='Add consensus sequences for subtypes to probe design', action='store_true')
parser.add_argument('-st', '--subtype_threshold', help='Insert a value for a minimum number of sequences of a given subtype it takes to justify generating and including a subtype consensus', default=0)
parser.add_argument('-ot', '--only_typed', help='Only sequences with a known geno- and subtype will be used in probe design', action='store_true')
parser.add_argument('-eg', '--exclude_genotype', nargs='*', help='list of space-seperated genotypes to be excluded for the purposes of probe design, but still tested against designed probes. If -cg is on, a consensus for these sequence is still be used to generate initial probes.')

args = parser.parse_args()
references = args.references_path
seq_similarity = float(args.similarity)
subtype_threshold = int(args.subtype_threshold)
threads = int(args.threads)

if args.exclude_genotype:
	exclude_genotypes = (' ').join(args.exclude_genotype)
else:
	exclude_genotypes = ''


# Here follows a lot of housekeeping. Using the path to the references, a variety of directories and new paths to intermediary files, result files, and the log file are derived.
directory = os.path.split(references)[0] # The directory where the references are
os.makedirs(os.path.join(os.path.split(references)[0], 'Probe_design_preliminary'), exist_ok=True) # Create the folder 'Probe_design_preliminary' in the directory of the references
Probe_design_preliminary = os.path.join(os.path.split(references)[0], 'Probe_design_preliminary') # Define the path to this new folder
os.makedirs(os.path.join(os.path.split(references)[0], 'Results'), exist_ok=True) # Create the folder 'Results' in the directory of the references
results_directory = os.path.join(os.path.split(references)[0], 'Results') # Define the path to this new folder
log_path = os.path.join(results_directory, os.path.split(references)[1].split('.')[0] + '_log.txt') # The path to the log file is defined (in the Results directory)
refs_wo_IUPAC_path = os.path.join(Probe_design_preliminary, 'refs_wo_IUPAC.fasta') # refs_wo_IUPAC_path will be the file that contains the references without IUPAC code
cleaned_refs_path = os.path.join(Probe_design_preliminary, 'cleaned_refs.fasta') # cleaned_refs_path will be the file that contains the references without IUPAC code and 3' UTR
os.chdir(os.path.dirname(subprocess.list2cmdline(sys.argv[0:1]))) # I change the current directory to that of the KTC_probe_designer_master.py
# Since I just changed current directory to that of the launched script, if sub-scripts are in the same directory their paths can simply be defined as:
Consensus_generator_script = './KTC_Consensus_generator_muscle.py'
Probe_design_script = './KTC_Probe_designer_v3.4.py'

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

logger('\t-- Script initiated with the following command: --\n' + subprocess.list2cmdline(sys.argv[:]) + '\n') # logs the command as it was launched
logger('\t-- Automated probe designer initialized --\n')

# =============================================================================
# Cleaning degenerate IUPAC code from references
# =============================================================================
# Probes cannot contain IUPAC code. Therefore, I break down any IUPAC code to one of its constituent bases at random.
logger('\t-- Cleaning IUPAC from reference sequences --\n')

# References are stored in a list as Biopython records (contains sequence, ID, alphabet...)
refs = []
for seq_record in SeqIO.parse(references, 'fasta'):
	if args.only_typed:
		if len(str(seq_record.id).split('_')[0]) > 1:
			refs.append(seq_record)
	else:
		refs.append(seq_record)

# IUPAC code is defined in a dictionary with a list to each key of IUPAC code of its constituent bases.
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


ref_len = len(refs) # ref_len is the length of the variable 'refs' - it is the number of references
refs_total_length = 0 # Holds the total number of bases in references
IUPAC_count = 0 # Holds the total number of IUPAC code found
refs_wo_IUPAC = {} # This dictionary will hold the sequences of the references without IUPAC code. The keys will be the IDs of the references.

# Each reference is looped over with each character in its sequence that matches a character that is a key in the 'IUPAC' dictionary being replaced with a random base from the list of the dictionary with the key of that character.
import random
for ref_count, ref in enumerate(refs):
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
	logger_no_clocker(ref.id + ' cleaned ' '(' + str(ref_count+1) + '/' + str(ref_len) + ' references)')

logger('Cleaning of IUPAC completed - ' + str(IUPAC_count) + ' of ' + str(refs_total_length) + ' characters replaced (' + str(round((IUPAC_count/refs_total_length)*100, 2)) + '%)')

#The resulting 'cleaned' sequences are written to the output directory in the fasta format.
with open(refs_wo_IUPAC_path, 'w') as out:
		for ref in refs_wo_IUPAC:
			out.write('>' + ref + '\n' + refs_wo_IUPAC[ref] + '\n')

logger_no_clocker('')

# =============================================================================
# Removing the poly(T) tract from 3' ends of all references
# =============================================================================
# We cannot design probes that contain the long string of 3' Ts - we would capture mRNA and sequences of low complexity

logger("\t-- Excising polyT tract from 3\' of references --\n")

# The following loop reads through each of the last 400 nucleotides of every reference and sees if it can satisfy either of two conditions:
# - it finds two T's in a row AND at least 13 of the following 20 nucleotides are T's
# - it finds ten T's in a row
# If either happens, the sequence without the polyT region is set to the sequence up until one of the above circumstances are satisfied
# If either does NOT happen. The sequence is considered as not having a polyT region and the sequence without the poly(T) region is simply set as the entirety of the sequence
window = 400
cleaned_refs = []


def cleaner(refe):
	T_count = 0
	sequence = str(refe.seq).replace('-', '').upper()
	PolyT_found = False
	for n in range(len(sequence[-window:])): # Essentially: Do the loop 'window' number of times
		position = len(sequence)-window+n # At n = 0, this is the window'th last position, at n = window, this is the last position
		base = sequence[position]
		if base.upper() == 'T':
			T_count += 1
			if position+1 != len(sequence): # We do not need to check what the next Ts are if this is the last position
				if sequence[position+1] == 'T' and sequence[position:position+20].count('T') >= 13: # If the next base is a T AND the next 20 bases hold 13 Ts
					sequence_no_PolyT = sequence[0:position] # The sequence without Ts is everything up until the current position.
					logger_no_clocker('Suspicion of polyT region found and excluded for reference ' + refe.id + ':\n5\'...' + sequence_no_PolyT[-150:] + ' ▼ ' + sequence[position:] + ' \'3\n')
					PolyT_found = True
					break
				elif T_count == 10: # If 10 Ts have been found in a row
					sequence_no_PolyT = sequence[0:position-9]
					logger_no_clocker('Suspicion of polyT region found and excluded for reference ' + refe.id + ':\n5\'...' + sequence_no_PolyT[-150:] + ' ▼ ' + sequence[position-9:] + ' \'3\n')
					PolyT_found = True
					break
		else:
			T_count = 0
	if PolyT_found == False: # If over the course of scanning for evidence of 3' poly(T) tract none was found:
		logger_no_clocker('No suspicion of polyT region for reference: ' + refe.id + '\n5\'...' + sequence[-300:] + ' \'3\n')
		sequence_no_PolyT = sequence # The sequence without poly(T) is the sequence in its entirety
	cleaned_refs.append('>' + refe.id + '\n' + sequence_no_PolyT + '\n') # Append a string in the fasta format of the reference without poly(T)

for refe in SeqIO.parse(refs_wo_IUPAC_path, 'fasta'):
	if args.only_typed and len(str(refe.id).split('_')[0]) == 0:
		logger_no_clocker('Reference: ' + refe.id + ' was excluded since no subtype was found.')
	else:
		cleaner(refe)

# The cleaned references are written to a new fasta file
with open(cleaned_refs_path, 'w') as out:
	for cleaned_ref in cleaned_refs:
		out.write(cleaned_ref)


logger_no_clocker('')


# =============================================================================
# Generating initial probes
# =============================================================================
# The probe designer script initialized at the end of this script requires a set of probes to get it started. These probes should cover well the diversity of the input references.

initial_probes = {}

# Generating initial probes based on the consensus of all sequences:
if args.consensus:
	logger('\t-- Generating a consensus of all sequences --\n')
	command = 'python ' + Consensus_generator_script + ' ' + cleaned_refs_path + ' -nh' + ' -l ' + log_path # This line writes the command for launching the consensus generator script on a per subtype basis with the -nh option
	subprocess.run(command, shell=True)

	for sequence in SeqIO.parse(cleaned_refs_path.split('.')[0] + '_consensus.fasta', 'fasta'):
		consensus_str = str(sequence.seq)

	initial_probes['consensus'] = []
	for l in textwrap.wrap(consensus_str, 120):
		initial_probes['consensus'].append(l)
	
	initial_probes['consensus'] = initial_probes['consensus'][:-1]
	initial_probes['consensus'].append(consensus_str[-120:])

# First we divide each input reference into a folder based on its subtype defined by the header of the fasta sequence prior to the first '_'
logger('\t-- Dividing references based on subtype --')

ref_dict = {} # Dictionary that will contain the references with each key being a subtype
for ref in SeqIO.parse(open(cleaned_refs_path, 'r'), 'fasta'):
	label = ref.id.split('_')[0]
	if label not in ref_dict:
		ref_dict[label] = []
	ref_dict[label].append(ref)

logger_no_clocker('')

# The 'KTC_Consensus_Generator_muscle.py' script is called for each. It generates a consensus sequence in a new file for each subtype (see script)
logger('\t-- Generating consensuses of subtypes in references --')

sub_con_paths = [] # After the consensuses have been generated, the paths to their fasta files are stored here
subtype_count = str(len(ref_dict))
current_subtype_count = 0
for label in ref_dict:
	current_subtype_count += 1
	logger_no_clocker('')
	logger('Generating consensus of subtype ' + label + ' (' + str(current_subtype_count) + '/' + subtype_count + ')')
	SeqIO.write(ref_dict[label], os.path.join(Probe_design_preliminary, label + '.fasta'), 'fasta')
	command = 'python ' + Consensus_generator_script + ' ' + os.path.join(Probe_design_preliminary, label + '.fasta') + ' -nh' + ' -l ' + log_path # This line writes the command for launching the consensus generator script on a per subtype basis with the -nh option
	subprocess.run(command, shell=True)
	sub_con_paths.append(os.path.join(Probe_design_preliminary, label + '_consensus.fasta'))

logger_no_clocker('')

# The consensus of each subtype is written to the same file and the consensus-generating script is called again to generate a consensus of consensuses
logger('\t-- Generating a consensus of the consensuses --\n')
subtype_consensuses = [] # Read and store the consensuses as Biopython records
for sub_con_path in sub_con_paths:
	for con in SeqIO.parse(open(sub_con_path, 'r'), 'fasta'):
		subtype_consensuses.append(con)

all_sub_con_path = os.path.join(Probe_design_preliminary, 'all_consensuses.fasta') # The path to the consensus of consensuses is defined
SeqIO.write(subtype_consensuses, all_sub_con_path, 'fasta') # All consensuses are written to that file
command = 'python ' + Consensus_generator_script + ' ' + all_sub_con_path + ' -nh' + ' -l ' + log_path# This line writes the command for launching the consensus generator on the consensuses script with the -nh option 
subprocess.run(command, shell=True)

# The sequence of the consensus of consensuses is loaded into the con_of_cons variable
con_of_cons_path = os.path.join(Probe_design_preliminary, 'all_consensuses_consensus.fasta')
for con in SeqIO.parse(con_of_cons_path, 'fasta'):
	con_of_cons = str(con.seq) # The con_of_cons variable is the sequence of the consensus of the consensuses

logger_no_clocker('')

# By now we should of course not have a polyT tract in the consensus of consensuses (all the sequences that went into making it had theirs removed)
# Regardless, we confirm that it does not have a polyT tract
logger('\t-- Checking Consensus of consensuses for 3\' polyT tract --\n')

window = 400
T_count = 0
PolyT_found = False
for n in range(len(con_of_cons[-window:])):
	position = len(con_of_cons)-window+n
	base = con_of_cons[position]
	if base == 'T':
		T_count += 1
		if T_count == 10:
			con_of_cons_noPolyT = con_of_cons[0:position-9]
			logger('Suspicion of polyT region found and excluded:\n5\'...' + con_of_cons_noPolyT[-150:] + ' ▼ ' + con_of_cons[position-9:] + ' \'3\n')
			PolyT_found = True
			break
	else:
		T_count = 0
if PolyT_found == False:
	logger('No suspicion of polyT region for consensus of consensuses:\n5\'...' + con_of_cons[-200:] + ' \'3\n')
	con_of_cons_noPolyT = con_of_cons

logger_no_clocker('')

# With the consensus of consensuses we then simply split it into words of 120. Since the last word is never 120 nt, we delete the last word and make a new one with the last 120 characters of the sequence.
logger('\t-- Creating probes from consensus of consensuses --\n')





if args.consensus_of_consensus:
	initial_probes['con_of_cons'] = []
	for l in textwrap.wrap(con_of_cons_noPolyT, 120):
		initial_probes['con_of_cons'].append(l)
	
	initial_probes['con_of_cons'] = initial_probes['con_of_cons'][:-1]
	initial_probes['con_of_cons'].append(con_of_cons_noPolyT[-120:])

logger_no_clocker('')

logger('\t-- Dividing subtype consensuses into genotypes --\n')

sub_con_dict = {} # Dictionary that will contain the consensus sequences of each subtype with each key being a genotype
for ref in SeqIO.parse(open(all_sub_con_path, 'r'), 'fasta'):
	label = ref.id[0]
	if label not in sub_con_dict:
		sub_con_dict[label] = []
	sub_con_dict[label].append(ref)

logger_no_clocker('')


# All the subtypes of each genotype is written to a file in preparation of consensus generation
logger('\t-- Generating consensuses of each genotype based on the consensus of its subtypes --')

genotype_con_paths = [] # After the consensuses have been generated, the paths to their fasta files are stored here
genotype_count = str(len(sub_con_dict))
current_genotype_count = 0
for label in sub_con_dict:
	current_genotype_count += 1
	logger_no_clocker('')
	logger('\t-- Generating consensus of genotype ' + label + ' subtype consensuses (' + str(current_genotype_count) + '/' + genotype_count + ') --')
	SeqIO.write(sub_con_dict[label], os.path.join(Probe_design_preliminary, label + '_SubT_cons.fasta'), 'fasta')
	command = 'python ' + Consensus_generator_script + ' ' + os.path.join(Probe_design_preliminary, label + '_SubT_cons.fasta') + ' -nh' + ' -l ' + log_path # This line writes the command for launching the consensus generator script on a per subtype basis with the -nh option
	subprocess.run(command, shell=True)
	genotype_con_paths.append(os.path.join(Probe_design_preliminary, label + '_SubT_cons_consensus.fasta'))


logger_no_clocker('')

if args.consensus_of_genotypes:
	for genotype_con_path in genotype_con_paths:
		genotype = os.path.split(genotype_con_path)[1].split('_')[0]
		initial_probes[genotype] = []
		for sequence in SeqIO.parse(genotype_con_path, 'fasta'):
			for l in textwrap.wrap(str(sequence.seq), 120):
				initial_probes[genotype].append(l)
			initial_probes[genotype] = initial_probes[genotype][:-1]
			initial_probes[genotype].append(str(sequence.seq)[-120:])



initial_probes_path = os.path.join(Probe_design_preliminary, 'initial_probes.fasta')
probe_count = 0
with open(initial_probes_path, 'w') as probe_out:
	for consensus in initial_probes:
		for p_no, p in enumerate(initial_probes[consensus]):
			probe_count += 1

			probe_out.write('>Probe_' + str(probe_count) + '_' + consensus + '\n' + p + '\n')


# =============================================================================
# Launching the probe designer script itself
# =============================================================================

# As a last step before the probe designer script is run, we take all the consensus sequences that we made and all the cleaned references to a new file to serve as templates for probe design
logger('\t-- Writing file with consensuses and references --\n')

cons_and_refs = []

if args.consensus_of_genotypes:
	for GTcon_path in genotype_con_paths:
		for GTcon in SeqIO.parse(GTcon_path, 'fasta'):
			cons_and_refs.append(GTcon)

if args.consensus_of_subtypes:
	for STcon in SeqIO.parse(all_sub_con_path, 'fasta'):
		cons_and_refs.append(STcon)

for ref in SeqIO.parse(cleaned_refs_path, 'fasta'):
	cons_and_refs.append(ref)

cons_and_refs_path = os.path.join(Probe_design_preliminary, os.path.split(references)[1].split('.')[0] + '.fasta')
SeqIO.write(cons_and_refs, cons_and_refs_path, 'fasta')

logger_no_clocker('')

# Finally, probe design is initiated
# We call upon the script using the initial probes we generated, the combined file with all the consensuses and cleaned references with the sequence similarity score and number of threads to be used
logger('\t-- Starting probe design --\n')

command = 'python ' + Probe_design_script + ' ' + initial_probes_path + ' ' + cons_and_refs_path + ' ' + results_directory + ' -s ' + str(seq_similarity) + ' -t ' + str(threads) + ' -eg ' + exclude_genotypes
logger(command)
subprocess.run(command, shell=True)

# =============================================================================
# Probe design and QC has finalized
# =============================================================================

logger('\t-- Probe Design finished --\n')