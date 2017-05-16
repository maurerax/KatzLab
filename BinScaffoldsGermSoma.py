


from Bio import SeqIO,Seq
from Bio.SeqUtils import GC
import os, sys, re

#------------------------------ Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
   
#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Runs Usearch on Split OrthoMCL Databases ----------------------###
###########################################################################################

def check_args(args):
	if len(args) == 0:
		print color.BOLD+'\nThis script is intended to '+color.GREEN+'HELP'+color.END+color.BOLD\
		+' you determine '+color.CYAN+'POTENTIALLY Germline Scaffolds '+color.END+color.BOLD\
		+' from a given Fasta File\n\n'
		print 'You will need to have run AUGUSTUS (or other gene prediction Software) on '\
		'the Fasta File with '+color.GREEN+'Genomic Scaffolds'+color.END+color.BOLD+' and '\
		'have BLASTed Transcripts AGAINST the Scaffolds File too!\n\n'
		print 'Predictions could be done AFTER or BEFORE training the prediction software'\
		' for your organism (if BEFORE, use Bacterial or "close" relative gene models)\n\n'
		print color.RED+'Example Usage:\n\n\t'+color.CYAN +'katzlab$ python BinScaffoldsGermSoma.py'\
		' ChilodonellaGermAssembly.1kb.fasta Chilo_WTA_NBU.GermAssembly.97ID.Ungapped.25wrd.tsv Chilo.GenePreds.Ecolik12.gtf\n\n'
		print color.ORANGE+color.BOLD+'\t\tQuestions/Comments? Email Xyrus (author) at maurerax@gmail.com\n\n'+color.END
		sys.exit()
	if len(args) != 3:
		print color.BOLD+color.RED+'\n\nCheck that you have included all the necessary command line inputs!\n\n'+color.END
		sys.exit()
	else:
		fasta = args[0]
		tsv = args[1]
		gtf = args[2]
		if gtf.lower().endswith('gtf') != True:
			print color.BOLD+'\n\nCheck that you truly have a '+color.CYAN+'GTF file'+color.END\
			+color.BOLD+' and '+color.RED+'NOT a GFF file'+color.END+color.BOLD+'.\n\n These'\
			"are not equivalent files and this script is designed for AUGUSTUS's "+color.CYAN\
			+'GTF files\n\n'+color.END
		
	return fasta,tsv,gtf
		
def set_up_dict(set_fasta):
	inFasta = [i for i in SeqIO.parse(set_fasta,'fasta') if len(i.seq) >= 5000]
	scaffold_dict = {}
	for i in inFasta:
		scaffold_dict.setdefault(i.description,[]).append(str(len(i.seq)))
		scaffold_dict.setdefault(i.description,[]).append('%.2f' % GC(i.seq))
	return	scaffold_dict	

def parse_TranscriptomeHits(intsv_file, scaffold_dict):
	inTSV = list(set(['\t'.join(i.split('\t')[:2]) for i in open(intsv_file).read().split('\n') if i != '']))
	core_TSV_Hits = [i.split('\t')[-1] for i in inTSV]
	for key in scaffold_dict.keys():
		scaffold_dict[key].append(str(core_TSV_Hits.count(key)))
	return scaffold_dict
	
def parse_gtf(gtf_file, scaffold_dict):
	inGTF = [i for i in open(gtf_file).read().split('\n') if i != '' and i.split('\t')[2] == 'gene']
	core_GTF_Hits = [i.split('\t')[0] for i in inGTF]
	for key in scaffold_dict.keys():
		scaffold_dict[key].append(str(core_GTF_Hits.count(key)))
		scaffold_dict[key].append(str(int(scaffold_dict[key][-1])*2000))
		scaffold_dict[key].append('%.2f' % (int(scaffold_dict[key][-1])*100/float(scaffold_dict[key][0])))
	return scaffold_dict

def categorize_data(scaffold_dict):
	for key in scaffold_dict.keys():
		if float(scaffold_dict[key][5]) > 20:
			scaffold_dict[key].append('Probably Non-Germline')
		else:
			if int(scaffold_dict[key][2]) > 0:
				scaffold_dict[key].append('Potentially Supported Germline')
			else:
				scaffold_dict[key].append('Potentially UnSupported Germline')
	scaffold_list = sorted((k+'\t'+'\t'.join(v) for k, v in scaffold_dict.items()),key=lambda x: -int(x.split('\t')[1]))
	LikelySupportedGerm = len([i for i in scaffold_list if i.split('\t')[-1] == 'Potentially Supported Germline'])
	LikelyUnSupportedGerm = len([i for i in scaffold_list if i.split('\t')[-1] == 'Potentially UnSupported Germline'])
	NonGerm = len(scaffold_list) - (LikelySupportedGerm + LikelyUnSupportedGerm)
	print '\nSupportedGerm == '+str(LikelySupportedGerm)+'\n'
	print 'UnSupportedGerm == '+str(LikelyUnSupportedGerm)+'\n'
	print 'NonGerm == '+str(NonGerm)+'\n'
	print 'TotalScaffolds == '+str(len(scaffold_list))+'\n\n'
	return scaffold_list
	
	
def write_out(summarized_data,outFileName,set_fasta):
	with open(outFileName,'w+') as w:
		w.write('Scaffold Name\tScaffold Length\tGC Content\tTranscript Hits\tPredicted ORFs\tTotal Predicted ORF Length\tPredicted % Coding\tType of Scaffold\n')
		for i in summarized_data:
			w.write(i+'\n')
	
	supported_germ_ID = [i.split('\t')[0] for i in summarized_data if i.split('\t')[-1] == 'Potentially Supported Germline']
			
	supported_germ = [i for i in SeqIO.parse(set_fasta,'fasta') if i.description in supported_germ_ID]
	
	with open(outFileName.replace('.tsv','.BSupGerm.fasta'),'w+') as x:
		for i in supported_germ:
			x.write('>'+i.description+'_BsupGerm\n'+str(i.seq)+'\n')

def main():
	args = sys.argv[1:]	
	fasta_file,tsv_file,gtf_file = check_args(args)
	
	outputName = fasta_file.split('.fas')[0]+'.CategorizedGermSoma.tsv'
		
	dict_pass1 = set_up_dict(fasta_file)
	dict_pass2 = parse_TranscriptomeHits(tsv_file, dict_pass1)
	dict_pass3 = parse_gtf(gtf_file, dict_pass2)
	Summary_list = categorize_data(dict_pass3)
	write_out(Summary_list,outputName, fasta_file)
	
	
main()


		