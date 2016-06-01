#!/usr/bin/env python

"""Script will produce dotplots using EMBOSS's dotmatcher. You will need to have .FASTA
files that contain 2 sequences that you want to compare (no more, no less). Then it 
will make the dot plots which should be decently named, and will note the direction of 
plotting (both fwd or one rev) in the title. Note that the plots only go in one direction.
so be careful when looking at them!

Also, note the 'cline's as you can adjust windowsize and threshold (%ID) there!"""



## To use this script simply call it as so:
##
##			katzlab$ python dotplot.py 
##
## If you haven't used chmod to make this an executable (and put in your path), then this 
## script needs to be in the same folder as all of your fasta files that you are comparing!
##
## Otherwise simply be in the directory and just call the script.
##
## Original intent is for this to be used towards the end of the MIC/MAC comparison pipeline!






from Bio import SeqIO
import os



def dotmatcher(f):
	## Just appends the sequences into a list so they can be separated prior to using
	## dotmatcher. 
	fasList = []
	infile = open(f,'r')
	for seq in SeqIO.parse(infile,'fasta'):
		fasList.append(seq.upper())
	
	## Separates your files into two distinct temp files, seq1 and seq2
	if len(fasList) > 0:
		for i in range(len(fasList)):
			for j in range(len(fasList)):
				if i != j:
					out = open('seq1.fasta','w')
					out2 = open('seq2.fasta','w')
					out.write('>' + fasList[i].id  + '\n' + str(fasList[i].seq).upper() + '\n')
					out2.write('>' + fasList[j].id  + '\n' + str(fasList[j].seq).upper()+ '\n')
					out.close()
					out2.close()
	## The actual command line call to dotmatcher. You can make adjustments here
					cline1 = 'dotmatcher -goutfile ' + f + '_dotplot_fwd -asequence=seq1.fasta -bsequence=seq2.fasta -windowsize 25 -threshold 70 -graph cps' 
					os.system(cline1)

	## Being lazy, so this section is the same as before, but does the reverse complement
	## for one of the sequences (just to be safe!)
	## There is a purpose for this duplicate chunk and that's to name the reversed sequence 
	## so it is easier to work with, without guessing. However coordinates won't be accurate
	## for the reversed sequence, so be careful with that! (They are actually descending)
	
	if len(fasList) > 0:
		for i in range(len(fasList)):
			for j in range(len(fasList)):
				if i != j:
					out = open('seq1.fasta','w')
					out2 = open('seq2.fasta','w')
	## Line below is where the reversal happens, also notes it in the name of the sequence too
					out.write('>' + fasList[i].id  + '_reverse\n' + str(fasList[i].seq.reverse_complement()).upper() + '\n')
					out2.write('>' + fasList[j].id  + '\n' + str(fasList[j].seq).upper()+ '\n')
					out.close()
					out2.close()
	## Here is the actual command again. Settings can be changed here, but should be identical
	# to the previous usage for better comparisons
					cline2 = 'dotmatcher -goutfile ' + f + '_dotplot_rev -asequence=seq1.fasta -bsequence=seq2.fasta -windowsize 25 -threshold 70 -graph cps' 
					os.system(cline2)

	## Just a cleanup step to move all the plots into a single location
def cleanup():
	os.system('mkdir DotPlots/')
	os.system('mv *.ps DotPlots/')

def main():
	for f in os.listdir(os.curdir):
		if f.endswith('.fasta'):
			dotmatcher(f)
#	cleanup()

main()


###
## Script was developed by Xyrus! If you have any questions, you can email me at:
## maurerax@gmail.com... although you do at your own peril. Fairly straight forward script
## and most issues are likely able to be solved with a quick google search.
##
## Remember that this requires the EMBOSS package as well as python and Biopython
##			