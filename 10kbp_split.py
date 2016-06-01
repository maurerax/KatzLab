#!/usr/bin/env python

"""This script can be used to split large assembled sequences (> 10kbp by default) into
10kbp (or smaller) chunks. The original intent was to split large Megabase sized assembled
sequences into BLASTable chunks to check for the density of bacterial vs eukaryotic proteins
to more effectively rule out large contigs from being eukaryotic"""


## To use this script simply call it as so:
##
##			katzlab$ python 10kbp_split.py filename 
##
## IF you haven't used chmod to make this an executable (and put in your path), THEN this 
## script needs to be in the same folder as all of your fasta files that you are comparing!
##
## Otherwise simply be in the directory and just call the script.
##
## Original intent is for this to be used towards the end of the MIC/MAC comparison pipeline!


from Bio import SeqIO
import sys

filename = sys.argv[1]


infile = SeqIO.parse(filename,'fasta')


for i in infile:
	count = 0
	newlist = []
	while count < len(i.seq):
		if (len(i.seq) - count) >= 10000:
			newlist.append('>'+str(i.description)+'_'+str(count)+'_'+str(count+10000)+'\n'+str(i.seq)[count:count+10000]+'\n')
			count += 10000
		else:
			newlist.append('>'+str(i.description)+'_'+str(count)+'_'+str(len(i.seq))+'\n'+str(i.seq)[count:count+10000]+'\n')
			count += 10000

	with open(str(i.description)+'_10kb_Splits.fasta','a') as x:
		for i in newlist:
			x.write(i)
