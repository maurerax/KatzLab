import sys, os
from itertools import product

#------------------------------ Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   PURPLE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script -------------------------------#

def prep_data(tsv_file, out_name):
	inTSV = []
	os.system('mkdir '+out_name+'_Germ_dir/')
	if open(tsv_file).read().count('\r') != 0:
		inTSV = [i for i in open(tsv_file).read().split('\r') if i != '']
	else:
		inTSV = [i for i in open(tsv_file).read().split('\n') if i !='']
	print color.BOLD + '\n\nThere are '+str(len(inTSV))+' Lines in the BLAST output' + color.END
	if len(inTSV) > 150000:
			print color.BOLD + '\n\nDue to the large number of lines, this can take a while...' + color.END
	return inTSV



def check_coverage(inTSV_list):
	print color.BOLD + '\n\nExtracting relevant information from BLAST output' + color.END
	coverage_dict = {}
	line_dict = {}
	for i in inTSV_list:
		coverage_dict.setdefault('\t'.join(i.split('\t')[:2]),[]).append(int(i.split('\t')[3]))
	for k,v in coverage_dict.items():
		if sum(coverage_dict[k]) >= 0.7*int(k.split('\t')[0].split('_')[3]):
			pass
		else:
			coverage_dict.pop(k,None)
	for i in inTSV_list:
		if '\t'.join(i.split('\t')[:2]) in coverage_dict.keys():
			line_dict.setdefault('\t'.join(i.split('\t')[:2]),[]).append('..'.join(i.split('\t')[6:10]))
	for k, v in line_dict.items():
		v.sort(key=lambda x: int(x.split('..')[0]))
	return line_dict


def identify_Singletons(post_cov_dict, pass_counter):
	singleton_dict = {}
	for k, v in post_cov_dict.items():
		if len(v) == 1:
			if pass_counter == 0:
				value_pair = [int(i) for j in v for i in j.split('..')[:2]]
				diff = abs(value_pair[0] - value_pair[1])
				if diff >= 0.75*int(k.split('\t')[0].split('_')[3]):
					singleton_dict[k] = post_cov_dict[k]
					post_cov_dict.pop(k, None)
			else:
				transcript_length = int(k.split('\t')[0].split('_')[3])
				hit_length = abs(int(v[0].split('..')[0])-int(v[0].split('..')[-1]))
				if hit_length >= 0.75* transcript_length:
					singleton_dict[k] = post_cov_dict[k]
				else:
					pass
				post_cov_dict.pop(k, None)
	return post_cov_dict, singleton_dict
	
	
def identify_repetitive(post_single_dict):
	print color.BOLD + '\n\nIdentifying '+ color.RED +'Repetitive ' + color.CYAN\
	+'Transcript-Genome Hits'+ color.END
	junk_pairs =[]
	for k, v in post_single_dict.items():
		repetitive = []
		odd_overlap = []
		for n in range(len(v)-1):
			x_range=set(range(min(int(s) for s in v[n].split('..')[:2]), max(int(s) for s in v[n].split('..')[:2])))
			y_range=range(min(int(s) for s in v[n+1].split('..')[:2]), max(int(s) for s in v[n+1].split('..')[:2]))
			germ3_MDS1 = set(range(min(int(s) for s in v[n].split('..')[-2:]), max(int(s) for s in v[n].split('..')[-2:])))
			germ5_MDS2 = range(min(int(s) for s in v[n+1].split('..')[-2:]), max(int(s) for s in v[n+1].split('..')[-2:]))
			germ_overlap = len(germ3_MDS1.intersection(germ5_MDS2))
			if len(x_range.intersection(y_range)) >= 0.7*len(x_range) or len(x_range.intersection(y_range)) >= 0.7*len(y_range):
				repetitive.append(k)
			else:
				if len(x_range.intersection(y_range)) >= 30:
					odd_overlap.append(k)
				if germ_overlap >= 10:
					odd_overlap.append(k)
				else:
					pass
		if len(repetitive) >= 0.3 * len(post_single_dict[k]):
			junk_pairs.append(k)
			#print k
			post_single_dict.pop(k,None)
		elif k in odd_overlap:
			junk_pairs.append(k)
			post_single_dict.pop(k,None)
		else:
			pass
	return post_single_dict, junk_pairs


def identify_intron_only(NoRepeat_dict, pointer_length_min):
	print color.BOLD + '\n\nIdentifying '+ color.RED +'Intronic ' + color.CYAN\
	+'Transcript-Genome Hits'+ color.END
	min_point = int(pointer_length_min)
	for k, v in NoRepeat_dict.items():
		pos_to_update = []
		introns = 0
		for n in range(len(v)-1):
			prime3_MDS1 = int(v[n].split('..')[1])
			prime5_MDS2 = int(v[n+1].split('..')[0])
			germ3_MDS1_pos = [int(s) for s in v[n].split('..')[-2:]]
			germ5_MDS2_pos = [int(s) for s in v[n+1].split('..')[-2:]]
			nearest_germ_dists = sorted(product(germ3_MDS1_pos, germ5_MDS2_pos), key=lambda t: abs(t[0]-t[1]))[0]
			germ3_MDS1_range = set(range(min(germ3_MDS1_pos), max(germ3_MDS1_pos)))
			germ5_MDS2_range = range(min(germ5_MDS2_pos), max(germ5_MDS2_pos))
			germ_overlap = len(germ3_MDS1_range.intersection(germ5_MDS2_range))
			germ_MDS1_direction = (germ3_MDS1_pos[0]-germ3_MDS1_pos[-1])/abs(germ3_MDS1_pos[0]-germ3_MDS1_pos[-1])
			germ_MDS2_direction = (germ5_MDS2_pos[0]-germ5_MDS2_pos[-1])/abs(germ5_MDS2_pos[0]-germ5_MDS2_pos[-1])
			if abs(prime3_MDS1 - prime5_MDS2) <= min_point and abs(nearest_germ_dists[0]-nearest_germ_dists[-1]) < 200:
				if germ_MDS1_direction == germ_MDS2_direction:
					if introns == 0:
						pos_to_update.append((n,))
						pos_to_update[-1] += (n+1,)
					else:
						pos_to_update[-1]+= (n,)
						pos_to_update[-1] += (n+1,)
					introns += 1
				else:
					introns = 0
			elif (abs(prime3_MDS1 - prime5_MDS2)+1) == germ_overlap:
				if introns == 0:
					pos_to_update.append((n,))
					pos_to_update[-1] += (n+1,)
				else:
					pos_to_update[-1]+= (n,)
					pos_to_update[-1] += (n+1,)
				introns += 1
			elif (abs(prime3_MDS1 - prime5_MDS2)) == germ_overlap:
				if introns == 0:
					pos_to_update.append((n,))
					pos_to_update[-1] += (n+1,)
				else:
					pos_to_update[-1]+= (n,)
					pos_to_update[-1] += (n+1,)
				introns += 1
			elif germ_overlap != 0 and (abs(prime3_MDS1 - prime5_MDS2)-1) == germ_overlap:
				if introns == 0:
					pos_to_update.append((n,))
					pos_to_update[-1] += (n+1,)
				else:
					pos_to_update[-1] += (n,)
					pos_to_update[-1] += (n+1,)
				introns += 1
			else:
				introns = 0
		for i in pos_to_update:
			new_soma = v[min(i)].split('..')[0]+'..'+v[max(i)].split('..')[1]
			new_germ = v[i[0]].split('..')[-2]+'..'+v[i[-1]].split('..')[-1]
			new_coords = new_soma+'..'+new_germ
			NoRepeat_dict[k].append(new_coords)	
		pos_ready_to_discard = list(set(sum(pos_to_update,())))
		for index in sorted(pos_ready_to_discard, reverse=True):
			del NoRepeat_dict[k][index]
	return NoRepeat_dict


def identify_IESs(MDS_dict):
	print color.BOLD + '\n\nIdentifying and Classifying '+ color.CYAN +'Transcript-Genome Hits '\
		+ color.END + color.BOLD + 'that contain '+color.RED +'IES Sequences\n\n' + color.END
	scrambled = []
	non_scrambled = []
	missing_info = []
	for k, v in MDS_dict.items():
		v.sort(key=lambda x: int(x.split('..')[0]))
		values = []
		for n in range(len(v)-1):
			germ3_MDS1 = set(range(min(int(s) for s in v[n].split('..')[-2:]), max(int(s) for s in v[n].split('..')[-2:])))
			germ5_MDS2 = range(min(int(s) for s in v[n+1].split('..')[-2:]), max(int(s) for s in v[n+1].split('..')[-2:]))
			germ_overlap = len(germ3_MDS1.intersection(germ5_MDS2))
			if int(v[n].split('..')[1]) - int(v[n+1].split('..')[0]) >= 1:
				if (int(v[n].split('..')[-2]) - int(v[n].split('..')[-1])) < 0:
					values.append('-')
				else:
					values.append('+')
				if (int(v[n+1].split('..')[-2]) - int(v[n+1].split('..')[-1])) < 0:
					values.append('-')
				else:
					values.append('+')
				if (int(v[n].split('..')[-1]) - int(v[n+1].split('..')[-2])) < 0:
					if germ_overlap == abs(int(v[n].split('..')[-1]) - int(v[n+1].split('..')[-2])):
						pass
					else:
						values.append('-')
				else:
					values.append('+')
			else:
				pass
		if len(set(values)) > 1:
			scrambled.append(k)
			MDS_dict[k].append('scrambled')
		elif len(set(values)) == 1:
			non_scrambled.append(k)
			MDS_dict[k].append('nonscrambled')
		else:
			MDS_dict[k].append('missing_mds')
	return scrambled, non_scrambled, missing_info, MDS_dict


def check_missing_MDS(post_MDS_ID_dict, non_scrambled_data):
		#### Write function to identify transcripts with MISSING MDS!
	missing_MDS = []
	for k, v in post_MDS_ID_dict.items():
		if k in non_scrambled_data:
			for n in range(len(v)-2):
				if (int(v[n].split('..')[1]) - int(v[n+1].split('..')[0])) <= -2:
					missing_MDS.append(k)
					post_MDS_ID_dict[k].pop(-1)
					post_MDS_ID_dict[k].append('missing_mds')
		post_MDS_ID_dict[k] = list(set(v))
		post_MDS_ID_dict[k].sort()
	non_scrambled = [i for i in non_scrambled_data if i not in missing_MDS]
	return non_scrambled, missing_MDS, post_MDS_ID_dict

def keep_StrongestEvidence(some_dict):
	junky = {}
	for k, v in some_dict.items():
#		print k
		indexes = []
		for n in range(len(v)-2):
			arr1 = tuple([int(i) for i in v[n].split('..')[-2:]])
			arr2 =tuple([int(i) for i in v[n+1].split('..')[-2:]])
			vals = sorted(product(arr1, arr2), key=lambda t: abs(t[0]-t[1]))[0]
			if abs(vals[0]-vals[1]) <= 20:
				indexes.append(n)
				indexes.append(n+1)
		if (len(v) - len(set(indexes))) <= 2:
			junky[k] = some_dict[k]
			some_dict.pop(k, None)
		else:
			indexes = list(set(indexes))
			for index in sorted(indexes, reverse=True):
				del some_dict[k][index]
			if len(some_dict[k]) <= 2:
				junky[k] = some_dict[k]
				some_dict.pop(k, None)
	return some_dict, junky
			
def write_out(Final_MDS_dict,IntronLess_Singles,Intron_Singles, out_name):
	Flattened_Data = []
	Flat_Singleton_Intronless = []
	Flat_Singleton_Introns = []
	
	for k, v in Final_MDS_dict.items():
		Flattened_Data.append([k+'\t'+'\t'.join(i.split('..'))+'\t'+v[-1] for i in v if 's' not in i])
	Flattened_Data = sorted((i for j in Flattened_Data for i in j),key=lambda x: (int(x.split('\t')[1].split('_')[1]),int(x.split('\t')[0].split('_')[1]),int(x.split('\t')[2])))
	for k, v in IntronLess_Singles.items():
		Flat_Singleton_Intronless.append([k+'\t'+'\t'.join(i.split('..'))+'\tNoIntron' for i in v if 's' not in i])
	for k, v in Intron_Singles.items():
		Flat_Singleton_Introns.append([k+'\t'+'\t'.join(i.split('..'))+'\tWithIntron' for i in v if 's' not in i])
	Singletons = Flat_Singleton_Introns + Flat_Singleton_Intronless
	Singletons = sorted((i for j in Singletons for i in j),key=lambda x: (int(x.split('\t')[1].split('_')[1]),int(x.split('\t')[0].split('_')[1]),int(x.split('\t')[2])))
	
	scram_count = len([k for k in Final_MDS_dict.keys() if 'scrambled' == Final_MDS_dict[k][-1]])
	nscram_count = len([k for k in Final_MDS_dict.keys() if 'nonscrambled' == Final_MDS_dict[k][-1]])
	mscram_count = len([k for k in Final_MDS_dict.keys() if 'missing_mds' == Final_MDS_dict[k][-1]])
	single_count = len(set([i.split('\t')[0] for i in Singletons]))
	print color.BOLD + 'There were '+str(scram_count)+' scrambled transcripts\n\n' 
	print color.BOLD +'There were '+str(nscram_count)+' non-scrambled transcripts in the data set\n\n'
	print 'There were '+str(mscram_count)+' transcripts with missing MDSs in the data set\n\n' 
	print 'There were '+str(single_count)+' transcripts lacking evidence for IESs in the data set\n\n'
	print 'There were '+str(scram_count+nscram_count+mscram_count+single_count)+' transcripts'\
	' mapped to Putative Germline Contigs\n\n'+color.END

	
	with open(out_name+'_Germ_dir/'+out_name+'.Scram_Data.tsv','w+') as w:
		for i in Flattened_Data:
			if '\tscrambled' in i:
				w.write(i+'\n')
	with open(out_name+'_Germ_dir/'+out_name+'.Missing_Data.tsv','w+') as w:
		for i in Flattened_Data:
			if '\tmissing_mds' in i:
				w.write(i+'\n')
	with open(out_name+'_Germ_dir/'+out_name+'.NonScram_Data.tsv','w+') as w:
		for i in Flattened_Data:
			if '\tnonscrambled' in i:
				w.write(i+'\n')
	with open(out_name+'_Germ_dir/'+out_name+'.SingleTons_NoIntron.tsv','w+') as w:
		for i in Singletons:
			w.write(i+'\n')
			
			
def main():

	if len(sys.argv) == 4:
		tsv_file = sys.argv[1]
		out_name = sys.argv[2]
		min_pointer = sys.argv[3]
	elif len(sys.argv) == 3:
		tsv_file = sys.argv[1]
		out_name = sys.argv[2]
		min_pointer = 2
	else:
		print color.BOLD+color.RED+'\nExample usage:\n\n\t'+color.CYAN+'katzlab$ python'\
		+' MDS_IES_Hunt2.py SomeProcessedTSVFile.tsv DesiredOutputDirectoryAndName\n\n'+color.END
		print color.ORANGE + color.BOLD + '\t\tQuestions/Comments? Email Xyrus (author) at maurerax@gmail.com\n\n'+color.END
		sys.exit()
		
	inTSV = prep_data(tsv_file, out_name)
	line_dict = check_coverage(inTSV)
	line_dict, IntronLess_Singles = identify_Singletons(line_dict, 0)
	line_dict, repetitive_seqs = identify_repetitive(line_dict)
	line_dict = identify_intron_only(line_dict, min_pointer)
	line_dict_MultiMDS, Intron_Singles = identify_Singletons(line_dict, 1)
	scrambled_keys, nonscrambled_keys, missing_info_keys, line_dict_MultiMDS_updated = identify_IESs(line_dict_MultiMDS)
	non_scrambled_keys_final, missing_MDS_keys, Final_MDS_dict = check_missing_MDS(line_dict_MultiMDS_updated, nonscrambled_keys)
#	best_evidence, junky_data = keep_StrongestEvidence(Final_MDS_dict)
	write_out(Final_MDS_dict,IntronLess_Singles,Intron_Singles, out_name)

main()
