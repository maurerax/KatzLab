from Bio import SeqIO
import sys
import os


### Preps the tsv file and your fasta file as well, splitting them so they 
### are easier to work with later (can find a way to work around splitting though if need be)

def prep_tsv_and_fasta(intsvfilename,BigFasta):
### Opens the BlastN report (outformat 6)
	infile = open(intsvfilename,'r').read().split('\n')
	infile.pop(-1)
	
### Grabs the names of the queries and their hits and only keeps those hits/query combos that are present more than once
	query_hit = [i.split('\t')[0]+'\t'+i.split('\t')[1] for i in infile]
	multi = list(set([i for i in query_hit if query_hit.count(i) > 1]))

### Writes out the lines of the tsv file into a new tsv file, smaller that has just a single query and its hit
	for j in infile:
		with open(j.split('|')[3]+'_'+j.split('|')[-2]+'_BlastNout.tsv','a') as w:
			for i in multi:
				if j.startswith(i):
					w.write(j+'\n')

### Takes the info from the previous chunk and grabs only those MAC contigs/transcripts that have hits in multiple spots in the MIC
	inseq = SeqIO.parse(BigFasta,'fasta')
	inseq = [i for i in inseq]
	for seq in inseq:
		for name in multi:
			if name.split('\t')[0] == seq.description:
				with open(seq.description.split('|')[-2]+'_Single.fas','a') as w:
					w.write('>'+seq.description+'\n'+str(seq.seq)+'\n')




#	os.system('mkdir PointerPrep')
#	for filename in os.listdir(os.curdir):
#		if '.tsv' in filename:
#			os.system('move '+filename+' PointerPrep')
#		if '.fas' in filename:
#			os.system('move '+filename+' PointerPrep')
	

def grab_pointer_sequences(intsvfilename, MACcontigfile):
### Defining all the variables and opening the ncessary files...
	inseq = SeqIO.parse(MACcontigfile,'fasta')
	inseq = [i for i in inseq]	
	seq_count = len(inseq)
	intsv = open(intsvfilename,'r').readlines()
	
	coords = []
	pointer = []
	pointer_seq = []
	count = 0

### Grab the raw coordinates of the IES-MDS junctions in the context of the MAC
	for line in intsv:
		coords.append(line.split('\t')[6])
		coords.append(line.split('\t')[7])
### Sorts the list as the coordinates may be out of order...		
	coords.sort()
### Identifies and grabs the region of the MAC BlastN hits that overlap (e.g the pointers) -- pointers under 10bp are identified	
	while count < len(coords)-1:
		if int(coords[count+1])-int(coords[count]) < 10:
			pointer.append(coords[count]+'_'+coords[count+1])
		count += 1
	
### Uses the pointer coordinates to isolate the actual pointer sequences
	for i in pointer:
		pointer_seq.append(str(inseq[0].seq)[(int(i.split('_')[0])-1):int(i.split('_')[1])])

### Opens a new file and writes the pointer sequences to that file (name is: '>MacContigName_PointerCoord:PointerStart_PointerEnd')		
	with open(MACcontigfile+'_Pointer_Sequences.txt','a') as w:
		count2 = 0
		for i in pointer_seq:
			if i != '':
				w.write('>'+str(inseq[0].description)+'_PointerCoord:'+pointer[count2]+'\n'+i+'\n')
			count2 += 1

def main():
	if len(sys.argv) == 3:
		intsvfilename = sys.argv[1]
		BigFasta = sys.argv[2]
		prep_tsv_and_fasta(intsvfilename,BigFasta)
		for filename in os.listdir(os.curdir):
			if filename.endswith('BlastNout.tsv'):
				intsvfilename = filename
				MACcontigfile = filename.split('_')[0]+'_Single.fas'
				grab_pointer_sequences(intsvfilename,MACcontigfile)
		
	else:
		for filename in os.listdir(os.curdir):
			if filename.endswith('BlastNout.tsv'):
				intsvfilename = filename
				MACcontigfile = filename.split('_')[0]+'_Single.fas'
				grab_pointer_sequences(intsvfilename,MACcontigfile)
		

main()