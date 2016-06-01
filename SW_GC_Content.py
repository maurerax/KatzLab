#!/usr/bin/env python

## This script uses biopython so ensure that it is installed! ##
## To run the script, the fasta file needs to be in the same directory
## Command is:   katzlab$ python slidingGCWindow.py
## Answer the info as you're asked!

from Bio import SeqIO
import os

print "\n"
filename = raw_input("Give the name of your .fasta file.  ")
print "\n"
window = int(raw_input('What is the window size you wish to use? (e.g. 10bp, 50bp, 100bp)  '))
print "\n"
steps = int(raw_input('What is the step size? (e.g. every 5 bp, every 40 bp?)  '))




##Function allows for user defined length of window and step size (effectively cuts the sequence to these lengths/steps##
def chunks(seq, win,step):
	seqlen = len(seq)
	for i in range(0,seqlen,step):
		j = seqlen if i+win>seqlen else i+win
		yield seq[i:j]
		if j==seqlen: break


def main():
	
	for seq_record in SeqIO.parse(filename,'fasta'):
		seq = seq_record.seq
		seq.upper()
		##Ensure that the sequence is in CAPS (just in case), then opens the output and calculates the GC content##
		##of the windows (separated by whatever step size) which can be graphed in excel... ##
		with open(seq_record.description+'_GC_Content_SlidingWindow_'+str(window)+'_Step_'+str(steps)+'.csv','a') as w:
			count = 1
			for subseq in chunks(seq,window,steps):
				GCcontent = (subseq.count('G') + subseq.count('C'))/ float(len(subseq))
				w.write(str(GCcontent)+','+str(count*steps)+'\n')
				count += 1
## Moves the files into a folder named the same as your fasta file!##
## Hunt in that folder for your sequences. Each sequence will have its own spreadsheet##
	os.system('mkdir '+filename.split('.')[0]+'_'+str(window)+'_SlidingWindow/')
	os.system('mv *GC_Content_SlidingWindow*.csv '+filename.split('.')[0]+'_'+str(window)+'_SlidingWindow/')
	
	print 'Finished! BUT Check the notes in the script if you are interested in determining the relative position of abnormal GC content regions!'
	
	
main()



##################################
## Notes on GC Content position ##
##################################

# To calculate the relative position of any data point produced in the dataset, you will 
# want to graph the data, or look at the row number (it's essentially the same). On graphs
# you should be able to identify the data points that you're interested in.
#
# In excel, you do this by hovering your cursor over the graph and it tells you the point 
# (e.g. point 1345) which corresponds to that row as well.
#
# To get the position in bp of that data point, you simply multiply the STEP size 
# (which is in the name of your spreadsheet!), NOT the window size (that's just how 
# wide your average is), by that row number to get the position.
#
# In this case let's say the point was 1234 and the step size was 40 and window size was 10.
# To get the position of that data point you simply multiply 1234 by 40:
# 1234 * 40bp = 49360bp
# The region of interest is around the 49360th base pair of the sequence we gave this script!

#Script by Xyrus! E-mail if you have questions! maurerax@gmail.com