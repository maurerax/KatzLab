#!/usr/bin/env python
"""
#####################################################################################################################
## This script is intended for use with a csv file that has coverage estimates per basepair.			   ##
## (Trying to remember) Those measurements can be done with bbmap which is part of the BBTools package		   ##
## you will have to check that usage. Once you produce a csv file with coverages (with Coverage_Estimate.csv as    ##
## the extension), this script will produce a graph of the sliding window average of the coverage for a given 	   ##
## contig/scaffold.												   ##
##														   ##
## If you do not have pylab installed as a python library, you will need to do so!				   ##
######################################################################################################################
"""

## USAGE:
## katzlab$ python SW_CovAvg.py windowsize stepsize --> runs on ALL files with the right extension in a folder ...
##
## katzlab$ python Sw_CovAvg.py 100 20 



import os
import pylab
import sys

if sys.argv < 3:
	print 'error, open script for usage info'
else:
	win = int(sys.argv[1])
	step = int(sys.argv[2])


def chunks(filename, win,step):
	infile = open(filename,'r').read().split('\n')
	cov_list = [int(i.split(',')[-1]) for i in infile if i != '']
	cov_list_len = len(cov_list)
	for i in range(0,cov_list_len,step):
		j = cov_list_len if i + win > cov_list_len else i + win
		yield cov_list[i:j]
		if j == cov_list_len: break

def rolling_averages(filename,win,step):
	x = []
	y = []
	with open(filename.split('Coverage')[0]+'RollCovAverages_'+str(win)+'Win_'+str(step)+'Step.csv','a') as w:
		count = 1
		for subseq in chunks(filename, win, step):
			averages = (sum(subseq)/float(len(subseq)))
			y.append(averages)
			w.write(str(averages)+','+str(count*int(step))+'\n')
			x.append(count*int(step))
			count += 1
	## Plot generation below with Pylab (pretty sweet!)
		pylab.plot(x,y,'r')
		pylab.savefig(filename.split('Coverage')[0]+'RollCovAverages_'+str(win)+'Win_'+str(step)+'Step.png')
		pylab.close()
	## Update folder name, if you want ... this was for testing originally
		os.system('mv *Step.csv RollingAverages_Fake/')
		os.system('mv *.png RollingAverages/Graphs_Fake/')
		
		
def cleanup():
	## Not a fully complete cleanup step yet, unsure if it seems wise to delete all files after plot
	## generation ... so many files.
	## Also rename folders to your liking, but commit those changes to all instances in the script.
	os.system('mkdir RollingAverages_Fake/')
	os.system('mkdir RollingAverages/Graphs_Fake')


def main():
	cleanup()
	for filename in os.listdir(os.curdir):
		if filename.endswith('Coverage_Estimate.csv'):

			print "\nGraphing contig "+str(filename.split('_Coverage')[0])
			rolling_averages(filename,win,step)		

main()

###########
## NOTES ##
###########
## Script written by Xyrus. Usage info and what the intent is can be found at the top! If you have any other questions
## feel free to email: maurerax@gmail.com
