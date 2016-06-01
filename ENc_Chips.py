#!/usr/bin/env python

print 'This calculates the codon bias of your input file - change file name in script. '
print 'for information about chips, see http://emboss.toulouse.inra.fr/cgi-bin/emboss/help/chips'


import os, re 
from Bio import SeqIO

for f in os.listdir(os.curdir):
	if re.search('.fas',f):
		input = f
		infile = open(input,'r')
		inseq = SeqIO.parse(infile,'fasta')
		for seq in inseq:
		    outfile = open('test.fas','w')
		    outfile.write('>' + re.sub('\|','_',seq.id) + '\n' + str(seq.seq) + '\n')
		    outfile.close()
		    os.system('chips -seqall test.fas -outfile ' + seq.id + '.chips')
		os.system('rm test.fas')
		outfile = open(input + '_chips.txt','w')
		for f in  os.listdir(os.curdir):
			if f.split('.')[-1] == 'chips':
				infile = open(f,'r')
				for line in infile:
					if re.search('Nc',line):
						NC = line.split()[-1] 
				outfile.write(f + ',' + NC + '\n')
		outfile.close()
		os.system('rm *chips')		
