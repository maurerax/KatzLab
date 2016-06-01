#!/usr/bin/env python



import os
import pylab
import sys

if sys.argv < 3:
	print 'error'
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
		pylab.plot(x,y,'r')
		pylab.savefig(filename.split('Coverage')[0]+'RollCovAverages_'+str(win)+'Win_'+str(step)+'Step.png')
#		pylab.close()
		os.system('mv *Step.csv RollingAverages_Fake/')
		os.system('mv *.png RollingAverages/Graphs_Fake/')
		
		
def cleanup():
	os.system('mkdir RollingAverages_Fake/')
	os.system('mkdir RollingAverages/Graphs_Fake')


def main():
	cleanup()
	for filename in os.listdir(os.curdir):
		if filename.endswith('Coverage_Estimate.csv'):

			print "\nGraphing contig "+str(filename.split('_Coverage')[0])
			rolling_averages(filename,win,step)		

main()