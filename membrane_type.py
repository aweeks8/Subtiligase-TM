#Written for python2.7
#Annotating TM types of proteins in N terminomics datasets
#Last updated by Amy on February 21, 2020.
#Dependencies: csv, sys, re, Biopython. 
#Run using command line arguments like this:
#python membrane_type.py Prospector Uniprot
#where Prospector is a csv output file from Prospector that has been filtered for only Abu peptides. Column 0 = observation number, Column 1 = Acc #, Column 2 = Uniprot ID, Column 3 = Num Unique, Column 4 = % Cov, Column 5 = Best Disc Score, Column 6 = Best Expect Val, Column 7 = m/z, Column 8 = z, Column 9 = ppm, Column 10 = Prev AA, Column 11 = DB peptide, Column 12 = Constant Modification, Column 13 = Variable Modification, Column 14 = Fraction, Column 15 = RT, Column 16 = Start, Column 17 = Expect, Column 18 = # in DB, Column 19 = Protein Length, Column 20 = Protein MW, Column 21 = Species, Column 22 = Protein name
#For this script, it's only important that Acc # is in Column 1, Start is in Column 16, and Protein Length is in Column 19, and also that the file is already filtered for Abu peptides. Can fix this at some point by identifying column indices based on name so the file requirements aren't so specific.
#Uniprot is the Uniprot flat file (.txt or .dat)
#Output is a csv file containing membrane protein type annotations.
#Example: python membrane_type.py HEK293T_272_nonquant_uniqueproteins.csv uniprot-20190809.txt

import csv
import sys
import re


alist, blist, clist, elist, flist, glist = [], [],[],[], [], []

with open(sys.argv[1], "rb") as fileA:
	reader = csv.reader(fileA)
	for row in reader:
		alist.append(row[2]) #Uniprot ID
		clist.append(row[16]) #Start
		elist.append(row[19]) #Protein Length
		flist.append(row[11]) #DB peptide
		glist.append(row[10]) #Prev AA

#make a matrix with all of the observed cleavage events and protein lengths

cleavages = []

for i in range(len(alist)):
	cleavages.append([alist[i], clist[i], elist[i], flist[i], glist[i]])


seq = set()

#make a list of all of the unique proteins in the dataset
matchlist = []
for row in alist:
	if row not in seq:
		blist.append(row)
		matchlist.append(''.join((' ', row)))
		seq.add( row )

matchlist2 = []

for row in matchlist:
	matchlist2.append(row.strip())

print(len(matchlist2))
print(len(blist))
print(matchlist2)

#erase contents of existing current.txt file

open('current.txt', 'wb').close()


#get all of the swissprot records for proteins in the dataset to avoid going through the entire uniprot file a bunch of times

with open(sys.argv[2], "rb") as fileA:
	counter1 = 0
	counter2 = 0
	for line in fileA:
		if "ID" in line[0:2] and any(item in line[0:] for item in matchlist):
			counter1 += 1
			with open('current.txt', 'a') as fileB:
				fileB.write(line)
		elif "//" in line[0:2] and counter1 == 1:
			counter2 +=1
			with open('current.txt', 'a') as fileB:
				fileB.write(line)
			counter1 = 0
			counter2 = 0
		elif counter1 == 1 and counter2 == 0:
			with open('current.txt', 'a') as fileB:
				fileB.write(line)

#open filtered uniprot file and and use biopython to get domain information

from Bio import SwissProt as sp
handle = open("current.txt")
features = [(record.entry_name, record.comments) for record in sp.parse(handle)]
domains = []
lengths = []
signalpep = []
propep = []
transmem = []
gpi = []
locations = []

print(len(features))

#terms to look for in features

domain = 'TOPO_DOM'
signal = 'SIGNAL'
propeptide = 'PROPEP'
transmembrane = 'TRANSMEM'
lipid = 'LIPID'
location = 'SUBCELLULAR LOCATION'

for j in range(len(features)):
	if any(item in features[j][0] for item in matchlist2):
		#print item
		for i in range(len(features[j][1])):
			if location in features[j][1][i]:
				locations.append((str(features[j][0]),str(features[j][1][i]).translate(None, '(),')))

#list of proteins and their subcellular locations
newlist = []

#counters to print to terminal 

type_I = 0
type_II = 0
type_III = 0
type_IV = 0
multi = 0
gpi = 0
others = 0
single = 0

for row in matchlist2:
	for i in range(len(features)):
		if str(row) in features[i][0]:
			for j in range(len(features[i][1])):
				if 'SUBCELLULAR LOCATION' in features[i][1][j] and 'type I ' in features[i][1][j]:
					newlist.append((row, 'Type I'))
					type_I += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'type II ' in features[i][1][j]:
					newlist.append((row, 'Type II'))
					type_II += 1
	
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'type III ' in features[i][1][j]:
					newlist.append((row, 'Type III'))
					type_III += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'type IV ' in features[i][1][j]:
					newlist.append((row, 'Type IV'))
					type_IV += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Multi-pass' in features[i][1][j]:
					newlist.append((row, 'Multi-pass'))
					multi += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'GPI-anchor' in features[i][1][j]:
					newlist.append((row, 'GPI-anchor'))
					gpi += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Single-pass ' in features[i][1][j]:
					newlist.append((row, 'Type IV'))
					single += 1
					break

print 'Total proteins: %d' % len(blist)
print 'Type I: %d' % type_I
print 'Type II: %d' % type_II
print 'Type III: %d' % type_III
print 'Type IV: %d' % type_IV
print 'Multi-pass: %d' % multi
print 'GPI-anchor: %d' % gpi
print 'Single-pass: %d' % single

with open("output.csv", "wb") as fileC:
		writer = csv.writer(fileC, lineterminator = '\n')
		writer.writerows(newlist)
