#!/usr/bin/env python2.7


##########################################################################################
## This script is intended to parse through a FASTA file of putative CILIATE germline	##
## sequences and with the output of a tab-delimited BLAST report (outfmt 6) that has	##
## been processing by the MDS_IES_Hunt.py script, will grab intervening IES sequences.	##
## 																						##
## Specifically will grab the IES sequences between 2 CONSECUTIVE MDSs... not from 		##
## pointer to pointer!																	##
##																						##
##	CURRENTLY (V1.0) WRITTEN FOR BLEPHARISMA'S PUTATIVE SCRAMBLED GERMLINE LOCI!!!		##
##																						##
##	Check that you have:																##
## 1. BLASTed a SOMATIC genome or transcriptome against the putative GERMLINE fasta		##
## 2. Characterized the BLAST hits as Scrambled vs Non-Scrambled with MDS_IES_Hunt.py	##
##																						##
##																						##
## 								COMMAND Example Below									##
##  python ParseIES.py Bleph.WTA.Germ.97ID.OUT.Scrambled.tsv Blepharisma_WGA.10k.fasta	##												
##																						##
## 			E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																						##
##########################################################################################

### Above field needs updating...

import collections, sys, os, numpy as np, scipy.stats as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

#----------------------------- Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script -------------------------------#

###########################################################################################
###--------------- Breaks a Sequence into Small User-Defined Sub-Seqs ------------------###
###########################################################################################

def chunks(seq, win,step):
	seqlen = len(seq)
	
	for i in range(0,seqlen,step):
		j = seqlen if i+win>seqlen else i+win
		yield seq[i:j]
		if j==seqlen: break


###########################################################################################
###--------------- Calculates GC Content for ALL Given Windows/Steps -------------------###
###########################################################################################
		
def get_GC_by_pos(seq_list, window, steps):

	print color.BOLD + '\n\nCalculating GC Content Across MDS-IES Boundaries'+color.END

	pos_dict = {}

	for n in range(-30,30+steps,steps):
		pos_dict.setdefault(n,[])

	for seq in seq_list:
		count = 0
		range_length = float(60+steps)/(steps)

		for subseq in chunks(seq,window,steps):
			count += 1
	
			pos = ((int(count-range_length)*steps)+30)
			GCcontent = GC(subseq)

			pos_dict[pos].append(GCcontent)
						

	return pos_dict


###########################################################################################
###----------- Merges GC Content Data to Single Value for a Given Window/Step ----------###
###########################################################################################

def summarize_GC_content(seq_list, steps, window):

	print color.BOLD + '\n\nSummarizing GC Content Information Across MDS-IES Boundaries'+color.END

	Summary_GC_list =[]

	boundary_gc_dict = get_GC_by_pos(seq_list, window, steps)

	for k, v in boundary_gc_dict.items():
	
		Avg_GC = (sum(v)/float(len(v)))
		lowerq = st.t.interval(0.95,len(v)-1,loc=np.mean(v), scale=st.sem(v))[0]
		upperq = st.t.interval(0.95,len(v)-1,loc=np.mean(v), scale=st.sem(v))[-1]

		Summary_GC_list.append(str(k)+'\t'+str(Avg_GC)+'\t'+str(lowerq)+'\t'+str(upperq))

		Summary_GC_list.sort(key=lambda x: int(x.split('\t')[0]))

	return Summary_GC_list


###########################################################################################
###--------- Extracts the Useful Info from TSV Blast Reports and Grabs Sequences -------###
###########################################################################################

def process_tsv_fasta_files(mic_name_dict, seq_dict):

	print color.BOLD + '\n\nParsing MDS-IES boundary positions from TSV and Fasta Files'+color.END

	fiveprimeMDS_seqs = []
	threeprimeMDS_seqs = []
	
	for k, v in mic_name_dict.items():
		for coords in v:
			if int(coords.split('::')[0]) < int(coords.split('::')[1]):
				cut_fiveprime = int(coords.split('::')[0])-30
				cut_fiveprime_interior = int(coords.split('::')[0])+30
				fiveprimeMDS_seqs.append(seq_dict[k][cut_fiveprime:cut_fiveprime_interior])
	
				cut_threeprime = int(coords.split('::')[1])-30
				cut_threeprime_interior = int(coords.split('::')[1])+30
				threeprimeMDS_seqs.append(seq_dict[k][cut_threeprime:cut_threeprime_interior])	
	
			if int(coords.split('::')[0]) > int(coords.split('::')[1]):
				cut_fiveprime = int(coords.split('::')[0])-30
				cut_fiveprime_interior = int(coords.split('::')[0])+30
				fiveprimeMDS_seqs.append(str(Seq(seq_dict[k][cut_fiveprime:cut_fiveprime_interior]).reverse_complement()))			
	
				cut_threeprime = int(coords.split('::')[1])-30
				cut_threeprime_interior = int(coords.split('::')[1])+30
				threeprimeMDS_seqs.append(str(Seq(seq_dict[k][cut_threeprime:cut_threeprime_interior]).reverse_complement()))
	
	fiveprimeMDS_seqs = list(set([i for i in fiveprimeMDS_seqs if i != '']))
	threeprimeMDS_seqs = list(set([i for i in threeprimeMDS_seqs if i != '']))
	
	return 	fiveprimeMDS_seqs, threeprimeMDS_seqs		
	
			
###########################################################################################
###------------------- Writing Method to Output Data to Spreadsheets -------------------###
###########################################################################################

def write_out(fivep_data_summary, threep_data_summary, tsv_file, window, steps):

	print color.BOLD + '\n\nWriting Data Out to New Spreadsheets\n\n'+color.END
	

	if os.path.isdir('MDS_IES_Boundary_GC/') != True:
		os.system('mkdir MDS_IES_Boundary_GC/')

	with open('MDS_IES_Boundary_GC/'+tsv_file.split('.tsv')[0]+'.5Prime_MDS_GC_Content_'+str(window)+'wndw_'+str(steps)+'step_V2.tsv','w+') as w:
		w.write('Position\tAverage GC\tLower_95CI\tUpper_95CI\n')
		for line in fivep_data_summary:
			w.write(line+'\n')	
			
	with open('MDS_IES_Boundary_GC/'+tsv_file.split('.tsv')[0]+'.3Prime_MDS_GC_Content_'+str(window)+'wndw_'+str(steps)+'step_V2.tsv','w+') as x:
		x.write('Position\tAverage GC\tLower 95CI\tUpper 95CI\tTotal MDSs\n')
		for line in threep_data_summary:
			x.write(line+'\n')

###########################################################################################
###--------------- Main Call (Checks Arguments and Runs Other Functions)  --------------###
###########################################################################################

def main():
	
	if len(sys.argv) != 3 and len(sys.argv) != 5:
		print color.BOLD + '\n\nIncorrect number of command-line arguments. Usage is:\n\n'\
		+ color.DARKCYAN + '\t\t\t katzlab$ python MyBLASToutputTSV.tsv MyGermlineFastaFile.fasta\n\n'\
		+ color.END
		sys.exit()
	else:
		tsv_file = sys.argv[1]
		fasta_file = sys.argv[2]
		if len(sys.argv) ==5:
			window = int(sys.argv[3])
			steps = int(sys.argv[4])
		else:
			window = 3
			steps = 2
	
	seq_dict = {}
	working_mic_names = {}
	
	inFasta = [i for i in SeqIO.parse(fasta_file,'fasta')]
	for seq_rec in inFasta:
		seq_dict[seq_rec.description] = str(seq_rec.seq)

	print color.BOLD + '\n\nOpening and Parsing through '+color.DARKCYAN+tsv_file+color.END\
	+ color.BOLD + ', depending on the file size, this can take a while...\n\n'

	intsv = [line for line in open(tsv_file,'r').read().split('\n') if line != '' and line.split('\t')[1] in seq_dict.keys()]

	if len(intsv) == 1 and intsv[0].count('\r') > 0:
		intsv = [line for line in open(tsv_file,'r').read().split('\r') if line != '' and line.split('\t')[1] in seq_dict.keys()]
	else:
		pass
		

	for line in intsv:
		working_mic_names.setdefault(line.split('\t')[1],[]).append('::'.join(line.split('\t')[8:10]))

#	total_MDS = str(len([pos for val in working_mic_names.values() for pos in val]))

	fiveprimeMDS_seqs, threeprimeMDS_seqs = process_tsv_fasta_files(working_mic_names, seq_dict)
	fivep_summary = summarize_GC_content(fiveprimeMDS_seqs, steps, window)
	threep_summary = summarize_GC_content(threeprimeMDS_seqs, steps, window)
	write_out(fivep_summary, threep_summary, tsv_file, window, steps)


main()

##########################################################################################
###-------------------------------------- NOTES ---------------------------------------###
##########################################################################################
#
# Currently script is V4.0 as of 03/20/2017
#
# Designed to estimate GC content in a given sliding window over a 60bp total interval
#
# This interval is not user-defined... but can be adjusted by updating the following functions:
#
#		process_tsv_fasta_files(); get_GC_by_pos()
#
# Questions or concerns? Email Xyrus (author) at maurerax@gmail.com
#
# Last Updated 03/07/2017
			



