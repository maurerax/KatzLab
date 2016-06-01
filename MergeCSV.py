#!/usr/bin/env python

import csv
from collections import OrderedDict
import os

rowname = raw_input('\nWhat is the row header?  (e.g. Sequence name, Taxon, etc...)  ')

filenames = []
for filename in os.listdir(os.curdir):
	if '.csv' in filename:
		filenames.append(filename)
file2 = tuple(filenames)

data = OrderedDict()
fieldnames = []

for filename in file2:
	with open(filename, 'rU') as fp:
		reader = csv.DictReader(fp)
		fieldnames.extend(reader.fieldnames)
		for row in reader:
			data.setdefault(row[str(rowname)], {}).update(row)

fieldnames = list(OrderedDict.fromkeys(fieldnames))
with open("Merged_Sheets.csv", "w") as wp:
	writer = csv.writer(wp)
	writer.writerow(fieldnames)
	for row in data.itervalues():
		writer.writerow([row.get(field, '') for field in fieldnames])
		

		
os.system('mkdir ArchivedCSV')
for filename in filenames:
	os.system('mv '+filename+' ArchivedCSV')


in_file = "Merged_Sheets.csv"  
out_file = "Merged_Sheets_fixed.csv"  
  
row_reader = csv.reader(open(in_file, "rb"))  
row_writer = csv.writer(open(out_file, "wb"))  
  
first_row = row_reader.next()  
row_writer.writerow(first_row)  
for row in row_reader:  
    new_row = [val if val else "0" for val in row] + (["0"] * (len(first_row) - len(row)))  
    row_writer.writerow(new_row)
    
os.system('rm Merged_Sheets.csv')