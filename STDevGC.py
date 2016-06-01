#!/usr/bin/env python

## This script is used to identify data points that exceed several different standard deviations
## from the average.
##
## Requires BioPython and is intended as a follow-up to slidingGCWindow.py
##
## command is:   katzlab$ python STDevGC.py filename.csv


from Bio import SeqIO

import numpy as np
import sys

filename = sys.argv[1]



inseqs = SeqIO.parse(filename,'fasta')
in1 = [i for i in inseqs]
GC_in = [str(i.description)+';'+str(((str(i.seq.upper()).count('C')+str(i.seq.upper()).count('G'))/float(len(str(i.seq)))))

infile = open(filename,'r').read().split('\n')
infile.pop(-1)

single = []
two = []
three = []



#######################################################################################
## Calculates all the min/max values needed and puts the value/step number in a list ##
#######################################################################################

def get_raw_value(filename):
	infileValues = [float(i.split(',')[1]) for i in infile]
	infileAverage = sum(infileValues)/float(len(infileValues))
	infileSTDEV = np.std(infileValues)
	SingleSDAbove = infileAverage+infileSTDEV
	SingleSDBelow = infileAverage-infileSTDEV
	TwoSDAbove = infileAverage+2*infileSTDEV
	TwoSDBelow = infileAverage-2*infileSTDEV
	ThreeSDAbove = infileAverage+3*infileSTDEV
	ThreeSDBelow = infileAverage-3*infileSTDEV
	
	for i in infile:
	
		if float(i.split(',')[0]) < ThreeSDBelow:
			three.append(i)
		if float(i.split(',')[0]) > ThreeSDAbove:
			three.append(i)
	
		if float(i.split(',')[0]) < TwoSDBelow:
			two.append(i)
		if float(i.split(',')[0]) > TwoSDAbove:
			two.append(i)
		
		if float(i.split(',')[0]) < SingleSDBelow and i not in two:
			single.append(i)
		if float(i.split(',')[0]) > SingleSDAbove and i not in two:
			single.append(i)
			
			
	print '\nThree SD: '+ str(ThreeSDAbove)+' --- '+str(ThreeSDBelow)+'\n'
	print str(len(three)) +' points beyond 3 STDs\n'
	print 'Two SD: '+ str(TwoSDAbove)+' --- '+str(TwoSDBelow)+'\n'
	print str(len(two)) +' points beyond 2 STDs\n'
	print 'One SD: '+ str(SingleSDAbove)+' --- '+str(SingleSDBelow)+'\n'
	print str(len(single)+len(two)) +' points beyond 1 STDs\n'

def get_From_Fasta(filename):
	infileValues = [float(i.split(';')[1]) for i in GC_in]
	infileAverage = sum(infileValues)/float(len(infileValues))
	infileSTDEV = np.std(infileValues)
	SingleSDAbove = infileAverage+infileSTDEV
	SingleSDBelow = infileAverage-infileSTDEV
	TwoSDAbove = infileAverage+2*infileSTDEV
	TwoSDBelow = infileAverage-2*infileSTDEV
	ThreeSDAbove = infileAverage+3*infileSTDEV
	ThreeSDBelow = infileAverage-3*infileSTDEV
	
	for i in infile:
	
		if float(i.split(';')[0]) < ThreeSDBelow:
			three.append(i)
		if float(i.split(';')[0]) > ThreeSDAbove:
			three.append(i)
	
		if float(i.split(';')[0]) < TwoSDBelow:
			two.append(i)
		if float(i.split(';')[0]) > TwoSDAbove:
			two.append(i)
		
		if float(i.split(';')[0]) < SingleSDBelow and i not in two:
			single.append(i)
		if float(i.split(';')[0]) > SingleSDAbove and i not in two:
			single.append(i)
			
			
	print '\nThree SD: '+ str(ThreeSDAbove)+' --- '+str(ThreeSDBelow)+'\n'
	print str(len(three)) +' points beyond 3 STDs\n'
	print 'Two SD: '+ str(TwoSDAbove)+' --- '+str(TwoSDBelow)+'\n'
	print str(len(two)) +' points beyond 2 STDs\n'
	print 'One SD: '+ str(SingleSDAbove)+' --- '+str(SingleSDBelow)+'\n'
	print str(len(single)+len(two)) +' points beyond 1 STDs\n'



######################################################################################
## Writes out to a file the value/step number that was some Standard Deviation away ##
######################################################################################

		
def outwrite(filename):
	
	if single != []:
		with open(filename.split('_GC_Content_Sliding')[0]+'_SingleSDPoints.csv','a') as w:
			for i in single:
				w.write(i+'\n')
	if two != []:
		with open(filename.split('_GC_Content_Sliding')[0]+'_TwoSDPoints.csv','a') as w:
			for i in two:
				w.write(i+'\n')
	if three != []:
		with open(filename.split('_GC_Content_Sliding')[0]+'_ThreeSDPoints.csv','a') as w:
			for i in three:
				w.write(i+'\n')

def main():
	#get_raw_value(filename)
	get_From_Fasta(filename)
	
	outwrite(filename)
	
	
main()
	
## Notes on this script!
#
# Intended use is after running slidingGCWindow.py as the input needs to be 
# organized as value, step size. The output is simply those rows in the spreadsheet that 
# satisfy whichever SD they fit into. The only redundancy allowed is for the 3rd Standard
# Deviation, which might be rare. Otherwise, each row in the csv file is assigned into its
# appropriate csv file.
#
# Script by Xyrus! E-mail if you have questions! maurerax@gmail.com
