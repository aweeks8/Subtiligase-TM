#Written for python2.7
#Script for adjusting pepXML output of Protein Prospector for import into Skyline
#Copyright (C) 2020 Amy M. Weeks

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#There is a discrepancy in the mass assigned to the N-terminal Abu modification by Protein Prospector and Skyline. This script modifies the Protein Prospector output so that it matches what is expected by Skyline.
#Run the script from the command line as follows:
#python adjustpepXML.py input output
#Example:
#python adjustpepXML.py Q20190719-09.pep.xml Q20190719-09_adj.pep.xml


import sys

f = open(sys.argv[1], 'r+b')
filedata = f.read()
f.close()
newdata= filedata.replace("""85.05276""", """86.15276""")
f = open(sys.argv[2], 'wb')
f.write(newdata)
f.close()



