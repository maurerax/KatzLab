#!/usr/bin/env python

## Script will calculate the global GC content of given sequences (rather than a sliding window)

from Bio import SeqIO
import sys

if len(sys.argv) < 4:
	queryfile = raw_input('What is your filename?  ')
	outname = raw_input('What do you want the output named?  ')
	Ngen = raw_input('Was the assembler NGen or SPAdes?   ')
else:
	queryfile = sys.argv[1]
	outname = sys.argv[2]
	Ngen = sys.argv[3]

if Ngen.lower() == 'spades':
	with open(outname+'.tsv','a') as w:
		qseqs = SeqIO.parse(queryfile, 'fasta')
		w.write('Sequence Name\tGC content\tLength\tCoverage\n')
		for seq_record in qseqs:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\t'+str(seq_record.description.split('_')[-3])+'\n')

if Ngen.lower() == 'ngen':

	with open(outname+'.tsv','a') as w:
		qseqs = SeqIO.parse(queryfile, 'fasta')
		w.write('Sequence Name\tGC content\tLength\tCoverage\n')
		for seq_record in qseqs:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\t'+str(int(seq_record.description.split('_')[1])*150/int(seq_record.description.split('_')[-2]))+'\n')
