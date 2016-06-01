#!/usr/bin/env python

"""
This script is intended for use with the output from read mapping with bbmap which is
part of the bbtools package by Brian Brushnell: https://sourceforge.net/projects/bbmap/   

########################################################################################
## OF NOTE:the isolate_names function can be used for any time you wish to separate a ##
##.csv file based on the name of the first column ()								  ##
########################################################################################
"""

## Must have run bbmap on your contigs prior to using this script (for now).
## Benchmark times:   2 contigs (~50kbp each) -- 100s; 12000 contigs (variable size) ~ 1200s
##
## To call this script with the raw output from BBMAP (see notes at end for bbmap command):
##
## katzlab$ SlidingWindow_CoverageEstimator.py TestBaseCov_AllTelos_MIC_Masurca convert
##
## If you already have a .csv and you wish to split it based on the first column:
##
## katzlab$ SlidingWindow_CoverageEstimator.py TestBaseCov_AllTelos_MIC_Masurca




import sys
import os
import re
import csv
from itertools import groupby
import pylab

#######################################################################################
## To make life easier (being safe), BBMAP outputs a tsv, and this function converts ##
## the file to a csv for the next function (not sure if it can handle tsv/I'm lazy)  ##
#######################################################################################


def convert_tsv(filename):
	infile = open(filename, 'r').readlines()
	with open(filename + '.csv','a') as w:
		for line in infile:
			line2 = re.sub('\t',',',line)
			w.write(line2)


#########################################################################################
## Opens the output of bbmap. Grabs the names of the contigs that were mapped and then ##
## separates the data so that each contig mapped has its own file with mapping data	   ##
#########################################################################################


def isolate_names(filename):
	for key, rows in groupby(csv.reader(open(filename)), lambda row: row[0]):
		print "\n Isolating " + key + "from " + str(filename.split('.')[0])
		with open("%s_Coverage_Estimate.csv" % key, 'w') as outcsv:
			for row in rows:
				outcsv.write(','.join(row) + '\n')
		os.system('mv *Coverage_Estimate.csv CoverageEstimates/')
	

############################################################
## Just sets up a folder for the estimates to be put into ##
############################################################

	
def cleanup1():
	os.system('mkdir CoverageEstimates/')


def main():
	if len(sys.argv) < 2:
		filename = raw_input('What is the name of your mapping stats file?  (e.g. TestBaseCov_AllTelos_MIC_Masurca)   ')
	if len(sys.argv) == 2:
		filename = sys.argv[1]
	if len(sys.argv) == 3:
		filename = sys.argv[1]
		convert_tsv(filename)
	cleanup1()
	isolate_names(filename+'.csv')



main()
	 



######################
## Notes on Script  ##
######################

# As mentioned before this script is designed for use AFTER mapping reads with BBMAP
# The intent is to look at per contig coverage disparities to identify drastic differences 
# in read coverage that may be biologically informative.
#
# Output is a .tsv file FOR EACH contig ... beware, it might be better to process contigs
# in smaller batches...
#
# Questions? Email Xyrus at maurerax@gmail.com

"""
# BBMAP commands (in order with examples):
#
# katzlab$ ./bbmap.sh ref=ChiloMIC_Masurca_5_8_2016_All_Telomeres.fasta
# 
# ChiloMIC_Masurca_5_8_2016_All_Telomeres.fasta is the file with the contigs
#
# katzlab$ ./bbmap.sh in=ChiloMac_Fwd103.fastq.gz basecov=TestBaseCov_AllTelos_MIC_Masurca
#
# ChiloMac_Fwd103.fastq.gz == reads to be mapped; TestBaseCov_AllTelos_MIC_Masurca == output with mapping data
"""