#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) < 5:
	queryfile = raw_input('What is your filename?  ')
	outname = raw_input('What do you want the output named?  ')
	length = raw_input('What is the minimum contig length?   ')
	Ngen = raw_input('Was the assembler NGen, SPAdes, Trinity or MaSuRCA?   ')
else:
	queryfile = sys.argv[1]
	outname = sys.argv[2]
	length = sys.argv[3]
	Ngen = sys.argv[4]


### Preps the fasta file for the subsequent steps AND writes out a fasta file 
### with the min contig length desired

InFasta = SeqIO.parse(queryfile,'fasta')
InFasta = InFasta
#PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
#print "\n\nThere are "+str(len(PrepFasta))+" contigs with size greater than "+str(int(length)-1)+"\n\n"
#out = open(outname+'.fasta','w')
#SeqIO.write(PrepFasta, out, 'fasta')



### Performs the stats that on the contigs that are above the minimum size that was selected
if Ngen.lower() == 'spades':
	PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	print "\n\nThere are "+str(len(PrepFasta))+" contigs with size greater than "+str(int(length)-1)+"\n\n"
	out = open(outname+'.fasta','w')
	SeqIO.write(PrepFasta, out, 'fasta')
	with open(outname+'.tsv','a') as w:
		w.write('Sequence Name\tGC content\tLength\tCoverage\n')
		for seq_record in PrepFasta:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\t'+str(seq_record.description.split('_')[-1])+'\n')

if Ngen.lower() == 'ngen':
	PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	print "\n\nThere are "+str(len(PrepFasta))+" contigs with size greater than "+str(int(length)-1)+"\n\n"
	out = open(outname+'.fasta','w')
	SeqIO.write(PrepFasta, out, 'fasta')
	with open(outname+'.tsv','a') as w:
		w.write('Sequence Name\tGC content\tLength\tCoverage\n')
		for seq_record in PrepFasta:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\t'+str(int(seq_record.description.split('_')[1])*150/int(seq_record.description.split('_')[-2]))+'\n')

if Ngen.lower() == 'masurca':
	InFasta = SeqIO.parse(queryfile,'fasta')
	InFasta = InFasta
	PrepFastaStats = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	PrepFastaOut = ['>'+seq.description+'_'+str(len(seq.seq))+'\n'+str(seq.seq)+'\n' for seq in PrepFastaStats]
	print "\n\nThere are "+str(len(PrepFastaStats))+" contigs with size greater than "+str(int(length)-1)+"\n\n"
	with open(outname+'.fasta','a') as x:
		for i in PrepFastaOut:
			x.write(i)
	with open(outname+'.tsv','a') as w:
		w.write('Sequence Name\tGC content\tLength\n')
		for seq_record in PrepFastaStats:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\n')

if Ngen.lower() == 'trinity':
	PrepFasta = [seq for seq in InFasta if int(len(seq.seq)) > (int(length)-1)]
	PrepFastaOut = ['>'+seq.description.replace(' ','_').split('_path')[0]+'\n'+str(seq.seq)+'\n' for seq in PrepFasta]
	print "\n\nThere are "+str(len(PrepFasta))+" contigs with size greater than "+str(int(length)-1)+"\n\n"
	out = open(outname+'.fasta','w')
	with open(outname+'.fasta','a') as x:
		for i in PrepFastaOut:
			x.write(i)
	with open(outname+'.tsv','a') as w:
		w.write('Sequence Name\tGC content\tLength\n')
		for seq_record in PrepFasta:
			w.write(str(seq_record.description)+'\t'+str(((str(seq_record.seq.upper()).count('C')+str(seq_record.seq.upper()).count('G'))/float(len(str(seq_record.seq)))))+'\t'+str(len(seq_record.seq))+'\n')
