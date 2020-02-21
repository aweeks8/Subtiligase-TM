#Written for python2.7
#Script for calculating distances between cleavage events and TM domains
#Copyright (C) 2020 Amy M. Weeks

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#Written by Amy on August 12, 2019. Last updated on February 21, 2020.
#Dependencies: matplotlib, numpy, biopython, csv, sys. 
#Run using command line arguments like this:
#python transmembranedistances.py Prospector Uniprot
#where Prospector is a csv output file from Prospector that has been filtered for only Abu peptides. Column 0 = observation number, Column 1 = Acc #, Column 2 = Uniprot ID, Column 3 = Num Unique, Column 4 = % Cov, Column 5 = Best Disc Score, Column 6 = Best Expect Val, Column 7 = m/z, Column 8 = z, Column 9 = ppm, Column 10 = Prev AA, Column 11 = DB peptide, Column 12 = Constant Modification, Column 13 = Variable Modification, Column 14 = Fraction, Column 15 = RT, Column 16 = Start, Column 17 = Expect, Column 18 = # in DB, Column 19 = Protein Length, Column 20 = Protein MW, Column 21 = Species, Column 22 = Protein name
#For this script, it's only important that Acc # is in Column 1, Start is in Column 16, and Protein Length is in Column 19, and also that the file is already filtered for Abu peptides. 
#Uniprot is the Uniprot flat file (.txt or .dat)
#Output is a csv file containing transmembrane distances and a pdf file containing a histogram of distances for cleavage events N-terminal to a transmembrane domain. Plotting parameters can be adjusted depending on the specific dataset as described below.
#Example: python transmembranedistances.py HEK293T_272_nonquant.csv uniprot-20190809.txt

#Imports

import csv
import sys
import re
import numpy as np
import matplotlib.pyplot as plt



alist, blist, clist, elist, flist, glist = [], [],[],[], [], []

#Read in Protein Prospector results file

with open(sys.argv[1], "rb") as fileA:
	reader = csv.reader(fileA)
	for row in reader:
		alist.append(row[2]) #Uniprot ID
		clist.append(row[16]) #Start
		elist.append(row[19]) #Protein Length
		flist.append(row[11]) #DB peptide
		glist.append(row[10]) #PrevAA

#make a matrix with all of the observed cleavage events and protein lengths

cleavages = []

for i in range(len(alist)):
	cleavages.append([alist[i], clist[i], elist[i], flist[i], glist[i]]) #Uniprot ID, Start, Protein Length, DB peptide, Prev AA


seq = set()

#make a list of all of the unique proteins in the dataset

for row in alist:
	if row not in seq:
		blist.append(row)
		seq.add( row )

#erase contents of existing current.txt file

open('current.txt', 'wb').close()

#get all of the swissprot records for proteins in the dataset to avoid going through the entire uniprot file a bunch of times


with open(sys.argv[2], "rb") as fileA:
	counter1 = 0
	counter2 = 0
	for line in fileA:
		if "ID" in line[0:2] and any(item in line[0:] for item in blist):
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
features = [(record.entry_name, record.features, record.sequence) for record in sp.parse(handle)]
location = [(record.entry_name, record.organelle) for record in sp.parse(handle)]
domains = []
lengths = []
signalpep = []
propep = []
transmem = []
gpi = []

domain = 'TOPO_DOM'
signal = 'SIGNAL'
propeptide = 'PROPEP'
transmembrane = 'TRANSMEM'
lipid = 'LIPID'

#features is a list of lists. features[j][i] contains the entry name, subsequents lists are derived from record.features created by Biopython.

for j in range(len(features)):
	if any(item in features[j][0] for item in blist):
		for i in range(len(features[j][1])):
			if domain in features[j][1][i]:
				domains.append((str(features[j][0]),str(features[j][1][i][1:2]).translate(None, '(),'), str(features[j][1][i][2:3]).translate(None, '(),'), str(features[j][1][i][3:4]).translate(None, '(),\'')))
			if signal in features[j][1][i]:
				signalpep.append((str(features[j][0]),str(features[j][1][i][1:2]).translate(None, '(),'), str(features[j][1][i][2:3]).translate(None, '(),'), str(features[j][1][i][3:4]).translate(None, '(),\'')))
			if propeptide in features[j][1][i]:
				propep.append((str(features[j][0]),str(features[j][1][i][1:2]).translate(None, '(),'), str(features[j][1][i][2:3]).translate(None, '(),')))
			if transmembrane in features[j][1][i]:
				transmem.append((str(features[j][0]),str(features[j][1][i][1:2]).translate(None, '(),'), str(features[j][1][i][2:3]).translate(None, '(),')))
			if lipid in features[j][1][i]:
				gpi.append((str(features[j][0]),str(features[j][1][i][1:2]).translate(None, '(),')))


#write csv file with protein names, domain boundaries, domain name, protein length. optional, can comment this out.

with open("currentoutput.csv", "wb") as fileC: 
		writer = csv.writer(fileC, lineterminator = '\n')
		writer.writerow(['name', 'start', 'end', 'domain', 'length', 'TM', 'TM start', 'TM end', 'distance in aa'])
		writer.writerows(domains)

#Determine if cleavage falls within domain boundaries, whether the protein has a TM domain, and how far away the cleavage is from the TM domain
for row in cleavages:
	for i in range(len(domains)):
		if row[0] in domains[i][0]:
			start = int(domains[i][1])
			end = int(domains[i][2])
			domainname = domains[i][3]
			cleavage = int(row[1])
			if cleavage >= start and cleavage <= (end + 1):
				row.append(domainname)
				row.append(start)
				row.append(end)
				for i in range(len(transmem)):
					if row[0] in transmem[i][0]:
						tmstart = int(transmem[i][1])
						tmend = int(transmem[i][2])
						dist = tmstart-int(row[1])
						row.append('Transmembrane')
						row.append(tmstart)
						row.append(tmend)
						row.append(dist)





#get position of transmembrane domain and determine where the cleavage event is relative to the transmembrane domain.

for row in cleavages:
	try:
		if row[5] == '':
			print 'empty'
		else:
			for i in range(len(transmem)):
				if row[0] in transmem[i][0]:
					start = int(transmem[i][1])
					end = int(transmem[i][2])
					cleavage = int(row[1])
					if cleavage >= start and cleavage <= (end + 1):
						row.append('Transmembrane')
						row.append(start)
						row.append(end)
						break
				else:
					row.append('-')
					break	
	except IndexError:
		pass

distances = []
seq = set()

#collect data for histogram. To be included, row must contain and Extracellular domain and a TM domain, and the cleavage must be on the N-terminal side of the TM domain.

for row in cleavages:
	try:
		if 'Extracellular' in row[5] and 'Transmembrane' in row[8] and row[11]>0 and row[12] == '-' and (str(row[0])+'_'+str(row[1])) not in seq:
			distances.append(row[11])
			unique = str(row[0])+'_'+str(row[1])
			seq.add( unique )
	except IndexError:
		pass


with open("tmdistances.csv", "wb") as fileD:
	writer = csv.writer(fileD, lineterminator = '\n')
	writer.writerows(cleavages)

#Plotting parameters. These can be adjusted depending on the specifics of the dataset.

histogram = plt.hist(distances, 200, facecolor = 'b')

plt.xlabel('TM distance')
plt.ylabel('number of cleavages')
plt.title('TM distance histogram')
plt.axis([0, 1000, 0, 50])
plt.savefig('histogram.pdf')


