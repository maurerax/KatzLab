#!/usr/bin/env python

## This script uses biopython so ensure that it is installed! 
## To run the script, the fasta file needs to be in the same directory
## Command is:   katzlab$ python SeqHunterSW.py
## Answer the info as you're asked!
## Also assumes that your hunted contigs came from SPAdes! Adjust as necessary!


##MAKE REVERSE COMPLEMENTS!!!! DUH! 



##############################################################################################
## Script will find the frequency that the windows generated from chunks(see the function)  ##
## occurs in the Hunter sequence. Then it writes these frequencies to a .csv file where row ##
## names refer to the 'hunter sequences' and their frequency in the hunted sequences...     ##
##############################################################################################


from Bio import SeqIO
import os

print "\n"
queryname = raw_input("Give the name of your .fasta file with sequences you are looking for (e.g. querying)  ")
print "\n"
subjname = raw_input('Give the name of your .fasta file that where query sequences will be hunted (e.g. subject)  ')
print '\n'
window = int(raw_input('What is the window size you wish to use? (e.g. 10bp, 50bp, 100bp)  '))
print "\n"
steps = int(raw_input('What is the step size? (e.g. every 5 bp, every 40 bp?)  '))


##Make empty dictionary for sequences names and their windows
xdict = {}
ydict={}
log = []




###########################################################
##Splits sequences to user defined window and step sizes ##
###########################################################
def chunks(seq, win,step):
	seqlen = len(seq)
	for i in range(0,seqlen,step):
		j = seqlen if i+win>seqlen else i+win
		yield seq[i:j]
		if j==seqlen: break
		
		
		
###################################################################
## Add windows to the dictionary that will be used for searching ##
###################################################################
def query_split(qname, win, step):
	for seq_record in SeqIO.parse(qname,'fasta'):
		seq = str(seq_record.seq.upper())
		key = seq_record.description
		xdict.setdefault(key,[])
		for subseq in chunks(seq, window, step):
			xdict[key].append(subseq)
			
###########################################################################################
## Takes reverse complement of seequences then uses them for searching (just like above) ##
###########################################################################################
def query_rcomp_split(qname, win, step):
	for seq_record in SeqIO.parse(qname,'fasta'):
		seqr = str(seq_record.seq.upper().reverse_complement())
		key = seq_record.description
		ydict.setdefault(key,[])
		for subseq in chunks(seqr, win, step):
			ydict[key].append(subseq) 
	
	
########################################################################################	
## Search through other .fasta file and calculates the frequency of the query windows ##		
########################################################################################
def sub_query_hunt(sname):

	for seq_record in SeqIO.parse(subjname,'fasta'):
		seq = str(seq_record.seq.upper())
		print str(seq_record.description)
		
		#### Formats and names CSV files based on SPAdes assembler output###
		with open(str(subjname)+'_'+str(window)+'Window_Log.txt','a') as w:
			for key, value in xdict.items():
				count = len([i for i in value if i in seq])
				if count != 0:
					w.write('F,'+str(key)+','+str(count/float(len(value)))+','+str(seq_record.description)+'\n')
					#log.append('F,'+str(key)+','+str(count/float(len(value)))+','+str(seq_record.description)+'\n')
					print str(str(key)+' : '+str(count/float(len(value)))+' : '+str(seq_record.description))
					
			for key, value in ydict.items():
				county = len([i for i in value if i in seq])
				if county != 0:
					w.write('R,'+str(key)+','+str(county/float(len(value)))+','+str(seq_record.description)+'\n')
					#log.append('R,'+str(key)+','+str(county/float(len(value)))+','+str(seq_record.description)+'\n')
					print 'R : '+str(str(key)+' : '+str(county/float(len(value)))+' : '+str(seq_record.description))
					
		#with open(seq_record.description.split('_cov')[0]+'_'+str(window)+'bp_SWFrequency_Step_'+str(steps)+'.csv','a') as w:
		#	w.write('Qseq,'+seq_record.description.split('_length')[0]+'\n')
		#	w.write(str(key)+','+str(count/float(len(value)))+'\n')

			
def main():
	query_split(queryname,window,steps)
	query_rcomp_split(queryname, window, steps)
	sub_query_hunt(subjname)
	#with open(str(subjname)+'HitsLog.txt','a') as w:
	#	for i in log:
	#			w.write(i)	
	print '\n\nFinished!\n'
	
	
	
main()

## Moves the files into a folder named the same as your fasta file!##
## Hunt in that folder for your sequences. Each sequence will have its own spreadsheet##
#os.system('mkdir '+subjname.split('.')[0]+'_Step_'+str(steps)+'_SeqHunterSW/')
#os.system('mv *Step_'+str(steps)+'.csv '+subjname.split('.')[0]+'_Step_'+str(steps)+'_SeqHunterSW/')



######################
## Notes on Script  ##
######################

# This script is designed to hunt for windows of a query (or many query sequences) in a 
# given .fasta file for searching.
#
# Ultimately designed for scanning the germline genome of ciliates for signatures of somatic
# genome sequences, particularly for alternatively processed gene loci (see Kovner and Katz 2010)
#
# Output is an excel spreadsheet for each sequence in the subject file (the one being hunted in)
#
#
# Script by Xyrus! E-mail if you have questions! maurerax@gmail.com