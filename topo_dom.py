#Written for python2.7

#Script for topological domain annotation of N terminomics datasets.
#Copyright (C) 2020 Amy M. Weeks

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#Last updated by Amy on February 21, 2020.
#Dependencies: matplotlib, numpy, biopython, csv, sys.#Run using command line arguments like this:
#python topo_dom.py Prospector Uniprot
#where Prospector is a csv output file from Prospector that has been filtered for only Abu peptides. Column 0 = observation number, Column 1 = Acc #, Column 2 = Uniprot ID, Column 3 = Num Unique, Column 4 = % Cov, Column 5 = Best Disc Score, Column 6 = Best Expect Val, Column 7 = m/z, Column 8 = z, Column 9 = ppm, Column 10 = Prev AA, Column 11 = DB peptide, Column 12 = Constant Modification, Column 13 = Variable Modification, Column 14 = Fraction, Column 15 = RT, Column 16 = Start, Column 17 = Expect, Column 18 = # in DB, Column 19 = Protein Length, Column 20 = Protein MW, Column 21 = Species, Column 22 = Protein name
#For this script, it's only important that Acc # is in Column 1, Start is in Column 16, and Protein Length is in Column 19, and also that the file is already filtered for Abu peptides. Can fix this at some point by identifying column indices based on name so the file requirements aren't so specific.
#Uniprot is the Uniprot flat file (.txt or .dat)
#Output is a csv file with topological domain annotations.
#Example: python topo_dom.py 20160206_Jurkat_untreated_degra_bothreps_uniquenterm_filtered.csv uniprot-human.txt

import csv
import sys
import re

#Read in the Protein Prospector results. 

alist, blist, clist, elist, flist, glist = [], [],[],[], [], []

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

matchlist = []

for row in blist:
	matchlist.append(''.join((' ', row)))

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

with open("currentoutput.csv", "wb") as fileC: #this is the section you were looking for.
		writer = csv.writer(fileC, lineterminator = '\n')
		writer.writerow(['name', 'start', 'end', 'domain', 'length'])
		writer.writerows(domains)

for row in cleavages:
	for i in range(len(domains)):
		if row[0] in domains[i][0]:
			start = int(domains[i][1])
			end = int(domains[i][2])
			domainname = domains[i][3]
			cleavage = int(row[1])
			if cleavage >= start and cleavage <= (end + 1):
				#print domainname
				row.append(domainname)
				row.append(start)
				row.append(end)
				break

for row in cleavages:
	try:
		if row[5] == '':
			print 'empty'
	except IndexError:
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

for row in cleavages:
	try:
		if row[5] == '':
			print 'empty'
	except IndexError:
		for i in range(len(signalpep)):
			try:
				if row[0] in signalpep[i][0]:
					start = int(signalpep[i][1])
					end = int(signalpep[i][2])
					cleavage = int(row[1])
					if cleavage >= start and cleavage <= (end + 1):
						row.append('Signal Peptide')
						row.append(start)
						row.append(end)
			except ValueError:
				pass

for row in cleavages:
	try:
		if row[5] == '':
			print 'hello'
	except IndexError:
		for i in range(len(propep)):
			if row[0] in propep[i][0]:
				try:
					start = int(propep[i][1])
					end = int(propep[i][2])
					cleavage = int(row[1])
					if cleavage >= start and cleavage <= (end + 1):
						row.append('Propeptide')
						row.append(start)
						row.append(end)
				except ValueError:
					pass

for row in cleavages:
	try:
		if row[5] == '':
			print 'hello'
	except IndexError:
		for i in range(len(signalpep)):
			try:
				if row[0] in signalpep[i][0]:
					start = int(signalpep[i][1])
					end = int(signalpep[i][2])
					cleavage = int(row[1])
					if cleavage >= start and cleavage <= (end + 1):
						row.append('Signal Peptide')
						row.append(start)
						row.append(end)
					else:
						row.append('Has signal peptide, but cleavage is not in it')
						for i in range(len(transmem)):
							try:
								if row[0] in transmem[i][0]:
									tmstart = int(transmem[i][1])
									if cleavage >= start and cleavage <= (tmstart + 1):
										row.append('Cleavage between signal peptide and TM domain')
										break
							except ValueError:
								pass
						for i in range(len(gpi)):
							try:
								if row[0] in gpi[i][0]:
									gpianchor = int(gpi[i][1])
									if cleavage >= start and cleavage <= (gpianchor + 1):
										row.append('Cleavage between signal peptide and GPI anchor')
										break
							except ValueError:
								pass
			except ValueError:
				pass

with open("output.csv", "wb") as fileD:
	writer = csv.writer(fileD, lineterminator = '\n')
	writer.writerow(['name', 'cleavage site', 'length', 'Abu peptide', 'P4-P1', 'Topo Domain', 'start', 'end'])
	writer.writerows(cleavages)