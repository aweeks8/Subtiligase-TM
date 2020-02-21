#Written for python2.7
#Script for annotation of subcellular locations of protein identified in N terminomics datasets
#Copyright (C) 2020 Amy M. Weeks

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#Written by Amy on December 11, 2017. Last updated on February 21, 2020.
#Dependencies: csv, sys, re
#Run using command line arguments like this:
#python test20171214.py Prospector Uniprot
#where Prospector is a csv output file from Prospector that has been filtered for only Abu peptides. Column 0 = observation number, Column 1 = Acc #, Column 2 = Uniprot ID, Column 3 = Num Unique, Column 4 = % Cov, Column 5 = Best Disc Score, Column 6 = Best Expect Val, Column 7 = m/z, Column 8 = z, Column 9 = ppm, Column 10 = Prev AA, Column 11 = DB peptide, Column 12 = Constant Modification, Column 13 = Variable Modification, Column 14 = Fraction, Column 15 = RT, Column 16 = Start, Column 17 = Expect, Column 18 = # in DB, Column 19 = Protein Length, Column 20 = Protein MW, Column 21 = Species, Column 22 = Protein name
#For this script, it's only important that Uniprot ID in the form XXX_HUMAN is in Column 2
#Uniprot is the Uniprot flat file (.txt or .dat)
#Output is two CSV files: one with subcellular locations, one with both subcellular locations and GO Cellular component information. The second file may contain multiple GO CC annotations per protein.
#Example: python subcellularlocation.py 20160206_Jurkat_lysate_untreated_bothreps.csv uniprot-human-2017.txt


#Generic/Built-in Imports

import csv
import sys
import re


alist, blist, clist, elist, flist, glist = [], [],[],[], [], []

with open(sys.argv[1], "rb") as fileA:
	reader = csv.reader(fileA)
	for row in reader:
		alist.append(row[2]) #Uniprot ID


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


#get all of the swissprot records for proteins in the dataset to avoid going through the entire uniprot file a bunch of times

open('current.txt', 'wb').close()

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

#open filtered uniprot file and and use biopython to get subcellular location information

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
human = '_HUMAN'

print len(features)
print features[1][0].index(human)


#terms to look for in features

domain = 'TOPO_DOM'
signal = 'SIGNAL'
propeptide = 'PROPEP'
transmembrane = 'TRANSMEM'
lipid = 'LIPID'
location = 'SUBCELLULAR LOCATION'
human = '_HUMAN'

for j in range(len(features)):
	for y in range(len(blist)):
		if features[j][0] == blist[y]:
			for i in range(len(features[j][1])):
				if location in features[j][1][i]:
					locations.append((str(features[j][0]),str(features[j][1][i]).translate(None, '(),')))


#list of proteins and their subcellular locations
newlist = []

#counters to print to terminal 

cell_membrane = 0
ER = 0
Mito = 0
Nuc = 0
Cyto = 0
Secreted = 0
membrane = 0
cell_surface = 0
others = 0

for y in range(len(blist)):
	for i in range(len(features)):
		if blist[y] == features[i][0]:
			for j in range(len(features[i][1])):
				if 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Cell membrane' in features[i][1][j]:
					newlist.append((blist[y], 'Cell membrane'))
					cell_membrane += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Membrane' in features[i][1][j]:
					newlist.append((blist[y], 'Membrane'))
					membrane += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Secreted' in features[i][1][j]:
					newlist.append((blist[y], 'Secreted'))
					Secreted += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Mitochondri' in features[i][1][j]:
					newlist.append((blist[y], 'Mitochondrion'))
					Mito += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Nucle' in features[i][1][j]:
					newlist.append((blist[y], 'Nucleus'))
					Nuc += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Cytoplasm' in features[i][1][j]:
					newlist.append((blist[y], 'Cytoplasm'))
					Cyto += 1
					break
				elif 'SUBCELLULAR LOCATION' in features[i][1][j] and 'Endoplasmic reticulum' in features[i][1][j]:
					newlist.append((blist[y], 'ER'))
					ER += 1
					break

print 'Total proteins: %d' % len(blist)
print 'Cell membrane: %d' % cell_membrane
print 'ER: %d' % ER
print 'Mitochondrion: %d' % Mito
print 'Nucleus: %d' % Nuc
print 'Cytoplasmic: %d' % Cyto
print 'Secreted: %d' % Secreted
print 'Membrane: %d' % membrane
print 'Others: %d' % others

with open("output_locations.csv", "wb") as fileC:
		writer = csv.writer(fileC, lineterminator = '\n')
		writer.writerows(newlist)

pmlist = []

for y in range(len(newlist)):
	with open('current.txt', "rb") as fileA:
		counter1 = 0
		counter2 = 0
		for line in fileA:
			if "ID" in line[0:2] and (''.join((' ', newlist[y][0]))) in line:
				protein = newlist[y][0]
				counter1 += 1
			elif "//" in line[0:2] and counter1 == 1:
				if counter2 == 0:
					pmlist.append((newlist[y][0], newlist[y][1]))	
					counter2 +=1
				counter1 = 0
				counter2 = 0
			elif counter1 == 1 and counter2 == 0:
				if "DR" in line[0:2] and "C:plasma membrane" in line:
					pmlist.append((newlist[y][0], newlist[y][1],'GO CC: plasma membrane'))
				elif "DR" in line[0:2] and "C:cell surface" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC: cell surface'))
				elif "DR" in line[0:2] and "C:integral component of plasma membrane" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC: integral component of plasma membrane'))
				elif "DR" in line[0:2] and "C:extracellular region" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:extracellular region'))
				elif "DR" in line[0:2] and "C:extracellular matrix" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:extracellular matrix'))
				elif "DR" in line[0:2] and "C:extracellular exosome" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:extracellular exosome'))
				elif "DR" in line[0:2] and "C:endoplasmic reticulum" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:endoplasmic reticulum'))
				elif "DR" in line[0:2] and "C:mitochondrion" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:mitochondrion'))
				elif "DR" in line[0:2] and "C:nucleus" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:nucleus'))
				elif "DR" in line[0:2] and "C:cytosol" in line:
					pmlist.append((newlist[y][0], newlist[y][1], 'GO CC:cytosol'))
				

with open("output_CC.csv", "wb") as fileC:
		writer = csv.writer(fileC, lineterminator = '\n')
		writer.writerows(pmlist)

