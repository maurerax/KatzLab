import sys
from Bio import SeqIO
import os

tsv_file = sys.argv[1]
query_fasta = sys.argv[2]
hit_fasta = sys.argv[3]




def cluster_info(filename):
	query_hit = []
	goodLines = []
	
	infile = open(filename,'r').read().split('\r')
			
	queries = [i.split('\t')[0] for i in infile if 'Telo' not in i]
	queries_multi = [i for i in queries if queries.count(i) > 1]
	hits = [i.split('\t')[2] for i in infile if 'Telo' not in i]
	hits_multi = [i for i in hits if hits.count(i) > 1]
	
	for i in queries_multi:
		for line in infile:
			if i == line.split('\t')[0]:
				goodLines.append(line)
				query_hit.append(i+':'+line.split('\t')[2])
			
	for j in hits_multi:
		for line in infile:
			if j == line.split('\t')[2]:
				goodLines.append(line)
				query_hit.append(line.split('\t')[0]+':'+j)

	goodLines = list(set(goodLines))
	query_hit = list(set(query_hit))

	return query_hit


def write_cluster_tsv(filename):
	with open(filename.split('tsv')[0]+'Multi_hits.tsv','a') as w:
		w.write(infile[0]+'\n')
		for i in goodLines:
			w.write(i+'\n')	

def bin_seqs(filename, qfile, hfile):		
	
	query_infile = SeqIO.parse(qfile,'fasta')
	hits_infile = SeqIO.parse(hfile,'fasta')
	
	query_infile = [i for i in query_infile]
	hits_infile = [i for i in hits_infile]
	
	for i in cluster_info(filename):
		with open(i.replace(':','_')+'_Bin.fasta','a') as w:
			for qseq in query_infile:
				if i.split(':')[0] == str(qseq.description.split('_')[0]):
					print qseq.description
					w.write('>'+qseq.description+'\n'+str(qseq.seq)+'\n')
			for hseq in hits_infile:
				if i.split(':')[1] in str(hseq.description):
					w.write('>'+hseq.description+'\n'+str(hseq.seq)+'\n')

def main():
	bin_seqs(tsv_file, query_fasta, hit_fasta)

main()

#	