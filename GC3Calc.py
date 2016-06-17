#!/usr/bin/env python

import string
from Bio import SeqIO
import sys


"""Usage has been updated. You must specify the file in the command line!"""
#  katzlab$ python GC3Calc.py myFastaFile.fasta


def checkorf(f,seq2id,nt):
	pos = {} #  eg pos[(43,T)] == 1 means the 43 nuc is a T and is in  codon position 1
	tripletdict = {} # triplet[54] = (char1,char2,char3) in codon 54
	firsts = []
	seconds = []
	thirds = []
	outfile = open(f.split('.')[0]+'_codonStats.txt','a')
	outfile.write(seq2id + '\t')
	flag = 0
	count = 0
	ccount = 0
	tripletdict[ccount] = []
	degensites = []
	for char in nt:
		count = count + 1				
		if flag == 0:
			pos[(count,char)] = 1
			tripletdict[ccount].append(char)	
			flag = 1
		elif flag == 1:
			pos[(count,char)] = 2
			tripletdict[ccount].append(char)	
			flag = 2
		elif flag == 2:
			pos[(count,char)] = 3
			tripletdict[ccount].append(char)	
			ccount = ccount + 1 
			tripletdict[ccount] = []
			flag = 0
	for double in pos.keys(): #double is position in seq, nucleotide in that position, eg (32,G) and 
		if pos[double] == 1:
			firsts.append(double[1])
		elif pos[double] == 2:
			seconds.append(double[1])
		elif pos[double] == 3:
			thirds.append(double[1])
			
	degensites = getGCdegen(tripletdict)		
	all = firsts + seconds + thirds
	print seq2id, str(len(all))
	#print seq2id + ':' + str(degensites.count('c') + degensites.count('C') + degensites.count('g') + degensites.count('G')) + ':' + str(len(degensites))
	
	GCall = float(all.count('c') + all.count('C') + all.count('g') + all.count('G'))/float(len(all))
	outfile.write('\tGCall\t' +  str(GCall))
	
	GCthirds = float(thirds.count('c') + thirds.count('C') + thirds.count('g') + thirds.count('G'))/float(len(thirds))
	outfile.write('\tGCthirds\t' +  str(GCthirds))
	
	
	GCdegen = float(degensites.count('c') + degensites.count('C') + degensites.count('g') + degensites.count('G'))/float(len(degensites))
	outfile.write('\tGCdegen\t' +  str(GCdegen) + '\tDegenSites\t'+str(len(degensites))+'\n')
	
	
	outfile.close()
	
def getGCdegen(tripletdict):
	degensites = []
	for codoncount in tripletdict.keys():
		codonstring = ''
		for char in tripletdict[codoncount]:
			if char in ['t','a','c','g','T','A','C','G']:
				codonstring = codonstring + string.upper(char)
			else:
				codonstring = codonstring + 'X'
		if len(codonstring) == 3:
			#print codonstring
			if codoncheck(codonstring) == '4-fold':
				degensites.append(tripletdict[codoncount][2]) #degensites = list of third positions at 4fold degenerate sites
	return degensites


def codoncheck(codon):
	ChiloCodonDict = {}
	ChiloCodonDict = {'TAA': 'Stop','TTA':'L','TTG':'L','CTT':'4-fold','CTC':'4-fold','CTA':'4-fold','CTG':'4-fold','GGT':'4-fold','GGC':'4-fold','GGA':'4-fold','GGG':'4-fold','TCT':'4-fold','TCC':'4-fold','TCA':'4-fold','TCG':'4-fold','AGT':'S','AGC':'S','ACT':'4-fold','ACC':'4-fold','ACA':'4-fold','ACG':'4-fold','CCT':'4-fold','CCC':'4-fold','CCA':'4-fold','CCG':'4-fold','CGT':'4-fold','CGC':'4-fold','CGA':'4-fold','CGG':'4-fold','AGA':'R','AGG':'R','GTT':'4-fold','GTC':'4-fold','GTA':'4-fold','GTG':'4-fold','GCT':'4-fold','GCC':'4-fold','GCA':'4-fold','GCG':'4-fold','TAT':'Y','TAC':'Y','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','ATT':'I','ATC':'I','ATA':'I','ATG':'M','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TTT':'F','TTC':'F'}
	if 'X' in codon:
		return 'bad codon'
	if ChiloCodonDict[codon] == '4-fold':
		return '4-fold'	
	if ChiloCodonDict[codon] == 'Stop':
		return 'Stop'

def main():
	if len(sys.argv) > 1:
		f = sys.argv[1]
	else:
		f = raw_input('What is your FASTA file named? ')
	inseq1 = SeqIO.parse(f,'fasta')
	for seq1 in inseq1:
		try:
			checkorf(f,seq1.id,str(seq1.seq))
		except:
			print seq1.id

main()
