#!/usr/bin/env python

## Script depends on BioPython
## This script will hunt for contigs/sequences in a .FASTA file that contain telomeres
## It will either grab all the sequences or just those with BOTH telomeres, given the input

from Bio import SeqIO
import os
import sys

######################################################################################
## Checks for the number of command line arguments. If not all of the arguments are ##
## present prompts the user to input them prior to running.							##
######################################################################################

if len(sys.argv) < 4:
	filename = raw_input('\nWhat is the name of the query fasta?  ')
	outfile = raw_input('What is the organism, assembler and date you assembled? e.g. ChiloMic_SPAdes_05_14_15  ')
	number = raw_input('Do you want ALL sequences with telomeres or those with BOTH?  ')
else:
	filename = sys.argv[1]
	outfile = sys.argv[2]
	number = sys.argv[3]

###########################################################################################
## Sets up the telomeric sequences (may adjust in future to accomodate other species...) ##
## as well as the lists that will capture the name of the contigs that have them.		 ## 
###########################################################################################

ref = ['CCCCAAACCCC','GGGGTTTGGGG']
id_list = []
id_revlist = []



###########################################################################################
## Generates two lists adding sequence names that have a telomere on either end. Kept as ##
## Separate lists in case user is interested in only one side of the sequences...        ##
###########################################################################################

def Telos(filename):
	for seq_record in SeqIO.parse(filename,'fasta'):
		for i in ref:
			if i in seq_record.seq[0:50]:
				id_list.append(seq_record.id)
			if i in seq_record.seq[-50:]:
				id_revlist.append(seq_record.id)
	


	
def main():
	seqiter = SeqIO.parse(filename, 'fasta')
	if number.upper() == 'ALL':
		Telos(filename)
		All = list(set(id_list+id_revlist))
		handle = open(outfile+'_AllTelomeres.fas','w')
		SeqIO.write((seq for seq in seqiter if seq.id in All),handle,'fasta')
		handle.close()
		print '\nThere are '+str(len(All))+' sequences with Telomeres\n'
	
	elif number.upper() == 'BOTH':
		Telos(filename)
		Both = [i for i in id_list if i in id_revlist]
		handle = open(outfile+'_BothTelomeres.fas','w')
		SeqIO.write((seq for seq in seqiter if seq.id in Both),handle,'fasta')
		handle.close()
		print '\nThere are '+str(len(Both))+' sequences with BOTH Telomeres\n'

	else:
		print '\n\nOne of your inputs does not make sense!\n\n Double check names and make sure you type either ALL or BOTH\n' 


main()

######################
## Notes on Script  ##
######################

# This script is designed solely for Chilodonella uncinata's telomere sequences
# If you want to use other telomeric repeats, make sure to change them in the 
# ref list, otherwise telomeric repeats are set already.
#
# Output will be a fasta file indicating the type of sequences that were searched for 
# (All vs. Both telomeres)
#
# Questions? Email Xyrus at maurerax@gmail.com