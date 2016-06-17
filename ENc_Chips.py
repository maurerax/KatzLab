#!/usr/bin/env python

print 'This calculates the codon bias of your input file - change file name in script. '
print 'for information about chips, see http://emboss.toulouse.inra.fr/cgi-bin/emboss/help/chips'


import os,re
from Bio import SeqIO
import sys

if len(sys.argv) > 1:
	infile = sys.argv[1]
else:
	infile = raw_input('What is your FASTA file named?  ')
inseq = SeqIO.parse(infile,'fasta')
for seq in inseq:
    outfile = open('test.fas','w')
    outfile.write('>' + re.sub('\|','_',seq.id) + '\n' + str(seq.seq) + '\n')
    outfile.close()
    os.system('chips -seqall test.fas -outfile ' + seq.id + '.chips')
os.system('rm test.fas')
outfile = open(infile + '_chips.txt','w')
for f in  os.listdir(os.curdir):
	if f.split('.')[-1] == 'chips':
		infile = open(f,'r')
		for line in infile:
			if re.search('Nc',line):
				NC = line.split()[-1] 
		outfile.write(f + ',' + NC + '\n')
outfile.close()
for filename in os.listdir(os.curdir):
	if filename.endswith('.chips'):
		os.system('rm '+filename)	
