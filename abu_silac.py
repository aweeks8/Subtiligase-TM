#Written for python2.7
#Script for SILAC data analysis.
#Copyright (C) 2020 Amy M. Weeks

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Written by Amy Weeks on April 12, 2019. Last updated February 21, 2020.
#Input files needed are: is a Skyline report containing Total Area MS1 values and idotp values for peptide in replicate datasets. See Supplementary Dataset 9 for formatting information AND a csv file containing all of the Protein Prospector reports for each experiment. See Supplementary Dataset 2 for example formatting.
#Output is a csv file with calculations appended to the Skyline report in subsequent columns.
#Run this script by entering the name of the Skyline report, the name of each replicate as in the Skyline report, and whether the experimental condition is the light or heavy isotope. Any number of experiment replicates can be accomodated.
#Example of how to run the script from the command line:
#python 20180412_abu_silac.py 20190803_652_HEK293T_278_Pervanadate_Abu_SILAC_quantification.csv Q20180422-04 light Q20180422-08 heavy Q20180422-10 heavy Q20180422-04_06_08_10.csv

import sys
import csv
import numpy as np
from scipy import stats
import math

alist = []

#Function to make a list of the datasets and experimental isotopes

def datasets(a,b):
	alist.append((a,b))

#Get filenames and experimental isotopes from sys.argv. In the future, may want to use a parameter file instead.

listlen = len(sys.argv)-1

for i in range(1,listlen):
	if i % 2 == 0:
		x = i
		y = i + 1
		z = sys.argv[y]
		j = sys.argv[x]
		datasets(j, z)
	else:
		pass

#Open Skyline SILAC quantitation file and append the results to a list.

abulist = []

with open(sys.argv[1], "rb") as fileA:
    reader = csv.reader(fileA)
    for row in reader:
        abulist.append(row)

#Find the indices of the columns that contain the light and heavy peak areas for each replicate; according the experimental isotope, calculate experiment to control ratio and append to row in abulist

ratiolist = []

for i in range(len(alist)):
	lightname = ''.join(('light', ' ', alist[i][0], ' ', 'Total Area MS1'))
	heavyname = ''.join(('heavy', ' ', alist[i][0], ' ', 'Total Area MS1'))
	lightindex = abulist[0].index(lightname)
	heavyindex = abulist[0].index(heavyname)
	rationame = (('%s log2 experiment-to-control ratio') % alist[i][0])
	abulist[0].append(rationame)
	ratioindex = abulist[0].index(rationame)
	ratiolist.append(ratioindex)
	for row in abulist[1:]:
		try:
			if alist[i][1] == 'light':
				exptctrlratio = math.log(float(row[lightindex])/float(row[heavyindex]), 2)
				row.append(exptctrlratio)
			elif alist[i][1] == 'heavy':
				exptctrlratio = math.log(float(row[heavyindex])/float(row[lightindex]), 2)
				row.append(exptctrlratio)
		except (ValueError, ZeroDivisionError):
			row.append('NA')
	
#Calculate mean ratio for each peptide. Null hypothesis is a log2 fold change of 0. Calculate pvalue using stats.ttest_ind. Calculate log2 ratios and -log10 pvalues.

abulist[0].append('log2 peptide mean ratio')
abulist[0].append('pvalue')
abulist[0].append('log10pvalue')
start = int(ratiolist[0])
end = int(start + len(ratiolist))
null = []
for i in range(len(ratiolist)):
	null.append(0)



for row in abulist[1:]:
	try:
		meanlist = []
		for i in range(start,end):
			meanlist.append(float(row[i]))
		mean = np.mean(meanlist)
		pvalue = stats.ttest_ind(meanlist, null, equal_var=True, nan_policy='propagate')
		row.append(mean)
		row.append(pvalue.pvalue)
		log10pvalue = -(math.log(float(pvalue.pvalue), 10))
		row.append(log10pvalue)
	except (ValueError, ZeroDivisionError):
		row.append('NA')
		row.append('NA')
		row.append('NA')

#gather information from Protein Prospector output and add it to the results report.

prospectorlist = []

with open(sys.argv[(len(sys.argv)-1)], "rb") as fileA:
	reader = csv.reader(fileA)
    	for row in reader:
        	prospectorlist.append(row)

abulist[0].append('Protein Preferred Name')
abulist[0].append('P4-P1')
abulist[0].append('start')
abulist[0].append('length')

for i in range(1, len(abulist)):
	for row in prospectorlist:
		if abulist[i][1] == row[11]:
			abulist[i].append(row[2])
			abulist[i].append(row[10])
			abulist[i].append(row[16])
			abulist[i].append(row[19])
			break

P4aa = abulist[0].index('P4-P1')
abulist[0].append('Color')


outputfile = sys.argv[1].replace('.csv', '_results.csv')

with open(outputfile, "wb") as fileD:
	writer = csv.writer(fileD, lineterminator = '\n')
	writer.writerows(abulist)
