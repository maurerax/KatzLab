#!/usr/bin/env python

##########################################################################################
## This script is intended to remove small transcripts or small contigs below a given	##
## minimum size from a transcriptome or genome assembly.								##
##																						##
## Prior to running this script, ensure the following:									##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 				##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder					##
##																						##
## 								COMMAND Example Below									##
##																						##
## 			E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																						##
##							Next Script(s) to Run: 										##
##	AutoBactVsEuk.py (removes SSU then Bact) or 2a_removeSSU.py then 2b_removeBact.py	##
##																						##
##########################################################################################


###--------------------------------------------------- Example Command ---------------------------------------------------###
#																															#
#	katzlab$ python ContigFilterPlusStats.py ../scaffolds.fasta LKH001_Spirostomum_rnaSPAdes_15_05_2016.300bp 300 spades	#
#																															#
#	Output data can be found in the 'LKH001_Spirostomum' folder; current scripts REQUIRE 'rnaSPAdes', update if needed		#					
#																															#
###-----------------------------------------------------------------------------------------------------------------------###


from Bio import SeqIO
import sys,os


#------------------------------ Colors For Print Statements ------------------------------#
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

#------------------------------ Checks the Input Arguments ------------------------------#

#if len(sys.argv) != 5 and len(sys.argv) !=4:
#	print color.BOLD + '\n\nOpen this script to double check that you have added all the necessary command-line inputs!\n\n'
#	print  color.RED + 'For example:\n\n\n' + color.BLUE + 'katzlab$ python ContigFilterPlusStats.py ../scaffolds.fasta LKH001_Spirostomum_rnaSPAdes_15_05_2016.300bp 300 spades\n\n' + color.END
#	sys.exit()

if len(sys.argv) == 5:
	queryfile = sys.argv[1]
	outname = sys.argv[2]
	length = sys.argv[3]
	Assembler = sys.argv[4]
	if 'bp' not in outname:
		print  color.BOLD + '\n\nYou must include the minimum length (with "bp") and "_rna" in the filename.\n\n\n'\
			+ color.RED + 'For example:\n\n' + color.BLUE + 'LKH001_Spirostomum_rnaSPAdes_15_05_2016.300bp\n' + color.END
		print color.BOLD + '\nWork-Flow is designed to recognize 300bp extentions by default (as a start)\n\n' + color.END
		sys.exit()
	else:
		pass
		
elif len(sys.argv) == 4:
	queryfile = sys.argv[1]
	outname = sys.argv[2]
	length = sys.argv[3]
	Assembler = 'none'
	if 'bp' not in outname:
		print  color.BOLD + '\n\nYou must include the minimum length (with "bp") and "_rna" in the filename.\n\n\n'\
			+ color.RED + 'For example:\n\n\n' + color.BLUE + 'LKH001_Spirostomum_rnaSPAdes_15_05_2016.300bp\n' + color.END
		print color.BOLD + '\nWork-Flow is designed to recognize 300bp extentions by default (as a start)\n\n' + color.END
		sys.exit()
	else:
		pass
	
#----------------------- Preparing Fasta File and Output Folders ------------------------#

InFasta = SeqIO.parse(queryfile,'fasta')
InFasta = [i for i in InFasta]

#----------- Renames the Contigs, Writes them out, and Calculates Basic Info -----------#

### Performs the stats that on the contigs that are above the minimum size that was selected
if Assembler.lower() == 'spades':
	
	PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	PrepFasta.sort(key=lambda seq: -len(seq.seq))
	
	print color.BOLD + "\n\nThere are " + color.BLUE + str(len(PrepFasta))+" contigs "\
		+ color.END + color.BOLD + " with size greater than " + color.RED + str(int(length)-1)\
		+"\n\n" + color.END

	with open(outname+'.fasta','w+') as x:
		count = 1
		for i in PrepFasta:
			x.write('>Contig_'+str(count)+'_length_'+str(len(i.seq))+'\n'+str(i.seq)+'\n')
			count += 1
			
	with open(outname+'.tsv','w+') as w:
		w.write('Sequence Name\tGC content\tLength\tCoverage\n')
		count = 1
		for seq_record in PrepFasta:
			w.write('Contig_'+str(count)+'_length_'+str(len(seq_record.seq))+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\t'+str(seq_record.description.split('_')[-3])+'\n')
else:
	count = 1

	PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	PrepFasta.sort(key=lambda seq: -len(seq.seq))

	print color.BOLD + "\n\nThere are " + color.BLUE + str(len(PrepFasta))+" contigs "\
		+ color.END + color.BOLD + " with size greater than " + color.RED + str(int(length)-1)\
		+"\n\n" + color.END	
	
	with open(outname+'.fasta','w+') as x:
		count = 1
		for i in PrepFasta:
			x.write('>Contig_'+str(count)+'_length_'+str(len(i.seq))+'\n'+str(i.seq)+'\n')
			count += 1
			
	with open(outname+'.tsv','w+') as w:
		count = 1
		w.write('Sequence Name\tGC content\tLength\n')
		for seq_record in PrepFasta:
			w.write('Contig_'+str(count)+'_length_'+str(len(seq_record.seq))+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\n')
			count += 1



#----------------------------------------- NOTES -----------------------------------------#
#
# This script is designed to take a given scaffolds file from the user and return JUST the
# assembled sequences greater than the given input as well as a summary spreadsheet with
# sequence lengths and GC content.
#
# If given an assembler name (e.g. SPAdes, MaSuRCA) it may provide more info in the output
# spreadsheet (like coverage)
#
# You will need:
#
# Biopython and an assembled transcriptome (or genome)
#
# This is the very FIRST! (or nearly last) script that is part of a bioinformatic
# pipeline designed for transcriptome assembly to transcriptome annotation/protein 
# prediction. 
#
# Example command:
#
# katzlab$ python 1_ContigFilterPlusStats.py scaffolds.fasta LKH001_Spirostomum_rnaSPAdes_15_05_2016.300bp 300 spades
#
# If you have any questions contact Xyrus (author): maurerax@gmail.com
#
# Updated 2/24/2017
