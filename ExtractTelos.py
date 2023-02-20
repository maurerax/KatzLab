#!/usr/bin/env python3.10

##__Updated__: 29_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python ExtractTelos.py
##__Options__: python ExtractTelos.py --help


from Bio import SeqIO
from Bio.Seq import Seq
import argparse, os, re, sys
from argparse import RawTextHelpFormatter,SUPPRESS



#------------------------------ Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   PURPLE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD+'\nThis script will extract '+color.ORANGE+'ALL'+color.END+color.BOLD+\
	' contigs from an assembly that are "capped" by at least one telomere.\n\nNote that '\
	'the telomeric sequence used needs to be provided!\n\nDefault outputs include fasta files'\
	' in a folder in the same directory of the Fasta File being processed.'+color.END+usage_msg(), 
	usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)

	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+" Fasta file of nucleotide sequences\n"+color.END)

	required_arg_group.add_argument('--telomere','-telo', default = 'CCCCAAT',
	help=color.BOLD+color.GREEN+" User-defined single Telomeric repeat\n (default = 'CCCCAAT')\n"+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	required_arg_group.add_argument('--minLen','-len', default=200, type=int,
	help=color.BOLD+color.GREEN+" Minimum contig length in base pairs\n (default = 200)"+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Print author contact information\n'+color.END)
	

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	quit_eval = return_more_info(args)
	if quit_eval > 0:
		sys.exit()

	args = parser.parse_args()
	
	args.outpath = os.path.abspath(args.input_file)
	
	out_folder_name = 'Telos_Checked_'+args.input_file.split('/')[-1].split('.fa')[0]
	
	args.out_folder = args.outpath.replace(args.input_file.split('/')[0],out_folder_name)
	
	args.both_seqs = args.input_file.split('/')[0].split('.fa')[0]+'.BothTelomeres.fasta'
	args.single_seqs = args.input_file.split('/')[0].split('.fa')[0]+'.SingleTelomere.fasta'
	
	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python ExtractTelos.py'\
	' --input_file Didinium_Genome_Assembly.fasta --telomere CCCCAAT '\
	'--minLen 200'+color.END


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		print (author)
		valid_arg += 1
	
	return valid_arg
	
###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	if os.path.isdir(args.out_folder) != True:
		os.system('mkdir '+args.out_folder)


###########################################################################################
###-------------------------- Preps Given Telomeric Sequence ---------------------------###
###########################################################################################

def prep_telo_seq(args):
	
	telo_seq = args.telomere
	half_pos_telo = int(round(len(telo_seq)/float(2)))
	telo_seq_prep = telo_seq*3

	fwd_telo_seqs = [telo_seq_prep[:len(telo_seq)+half_pos_telo], telo_seq_prep[half_pos_telo:len(telo_seq)+2*half_pos_telo]]
	rev_telo_seqs = [str(Seq(i).reverse_complement()) for i in fwd_telo_seqs]

	return fwd_telo_seqs, rev_telo_seqs


###########################################################################################
###----------------------- Grabs Contigs with Telomeric Sequence -----------------------###
###########################################################################################

def grab_telos(args):

	fwd_telos, rev_telos = prep_telo_seq(args)
	
	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]
	
	with_fwd = []
	with_rev = []

	print (color.BOLD+'\nHunting for contigs capped with telomeres using '+color.RED+\
	args.telomere+color.END+color.BOLD+' as the "seed"\n'+color.END)
	
	for i in inFasta:
		if len(re.findall(fwd_telos[0],str(i.seq)[:50])) > 0:
			with_fwd.append(i)
		if len(re.findall(fwd_telos[1],str(i.seq)[:50])) > 0:
			with_fwd.append(i)
		if len(re.findall(rev_telos[0],str(i.seq)[-50:])) > 0:
			with_rev.append(i)
		if len(re.findall(rev_telos[1],str(i.seq)[-50:])) > 0:
			with_rev.append(i)

# Removing duplicates

	# with_rev
	
	with_rev_single = []
	all_id_r = []

	for i in with_rev:
		if i.id not in all_id_r:
			with_rev_single.append(i)
			all_id_r.append(i.id)
				
	# with_fwd

	with_fwd_single = []
	all_id_f = []

	for i in with_fwd:
		if i.id not in all_id_f:
			with_fwd_single.append(i)
			all_id_f.append(i.id)

# Separating sequences with telomere markers at both ends from those with markers at only one end

	with_both = []
	single_telo = []

	for i in with_rev_single:
		if i.id in all_id_f:
			with_both.append(i)
		else:
			single_telo.append(i)


	print (color.BOLD+'\nOf the initial '+color.GREEN+str(len(inFasta))+color.END+color.BOLD+\
	' contigs, there were '+color.ORANGE+str(len(with_both))+color.END+color.BOLD+' complete'\
	' chromosomes (with both telomeres) and '+color.CYAN+str(len(single_telo))+color.END+\
	color.BOLD+' partial chromosomes (with a single telomere)\n\n'+color.END)
	
	return with_both, single_telo


###########################################################################################
###---------------------- Writes Data out to Separate Fasta Files ----------------------###
###########################################################################################
	
def save_telomeric_seqs(args):

	complete_chrom, partial_chrom = grab_telos(args)
	
	with open(args.out_folder+'/'+args.both_seqs, 'w+') as w:
		for i in complete_chrom:
			w.write('>'+i.description+'\n'+str(i.seq)+'\n')
	
	with open(args.out_folder+'/'+args.single_seqs, 'w+') as x:
		for i in partial_chrom:
			x.write('>'+i.description+'\n'+str(i.seq)+'\n')


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	args = check_args()
	
	prep_folders(args)
	
	save_telomeric_seqs(args)
	
main()
