"""
#Script for extracting sequences from BLAST database PLUS version 2, Used to compare ABD geneome of wheat using survey sequence
#Run this where blast database is installed (ie. location of [db.nhr, db.nin and db.nsq])
"""
#OUTPUT: report_x.csv and  report_xx.csv

import os
import csv
from subprocess import call
import sys

csvlist = []
dblist = []
buff = 500

#*************************PARSING OPTIONAL COMMAND LINE ARG***************
if len(sys.argv) > 1:
	try:
		buff = int(sys.argv[1])
	except (ValueError):
		print "USAGE! python " + sys.argv[0] + " <buffer size>"
		exit(1)

#*****************GATHERING AVAILABLE CSV FILES************************
for item in os.listdir("."):
    if item.endswith(".csv"):
        csvlist.append(item)

#*****************CSV SELECTION AND ERROR HANDLING*****************
print "Which csv file?"
for i in range(len(csvlist)):
	print '\t' + str(i + 1) + ": " + csvlist[i]

csvfile = raw_input("Enter number of selection: ")
filename = ""
try:
	csvfile = int(csvfile)
	if csvfile < 1:
		raise IndexError
	filename = csvlist[csvfile - 1]
except (ValueError, IndexError):
	print "Invalid number input"
	exit(1)

if not os.path.isfile(filename):
	print "Invalid file name"
	exit(1)

#*********************GATHERING AVAILABLE DABASES****************************
for item in os.listdir("."):
    if item.endswith(".nin") or item.endswith(".nsq") or item.endswith(".nhr"):
        dblist.append(item.replace('.nsq', '').replace('.nin', '').replace('.nhr', ''))

#**********************DATABASE SELECTION AND ERROR HANDLING******************
dblist = set(dblist)
dblist = list(dblist)

print "Which database?"
for i in range(len(dblist)):
	print '\t' + str(i + 1) + ": " + dblist[i]

dbnum = raw_input("Enter number of selection: ")
database = ""
try:
	dbnum = int(dbnum)
	if dbnum < 1:
		raise IndexError
	database = dblist[dbnum - 1]
except (ValueError, IndexError):
	print "Invalid number input"
	exit(1)

if '.' in database:
	database = database[:database.index('.')]

#*******PRINT BREAK
print 'Extracting FASTA...'

#************************SEQUENCE ADJUST METHOD********************
mapping = {'A/G':'R','G/A':'R', 'C/T':'Y','T/C':'Y', 'G/C':'S','C/G':'S', 'A/T':'W','T/A':'W', 'G/T':'K','T/G':'K', 'A/C':'M','C/A':'M'}
def seqAdjust(sequence, change, loc):
	tempseq = sequence[:loc-1] + mapping[change] + sequence[loc:]
	return tempseq

#***************************FASTAFY PART*************************************
textfile = open(filename, "rb")

outfile = 'tempfasta.fa'
f = open(outfile, 'w')

csvreader = csv.reader(textfile, delimiter=',', quotechar='"')
first = True
chromids = {} #For replacing BLORDS later
for line in csvreader:
	if first:
		first = False
	else:
		adjseq = seqAdjust(line[1].strip(), line[3].strip(), int(line[4].strip()))
		f.write(">" + line[0].strip() + "\n" + adjseq + "\n")
		#BLORD replacement prep
		tempkey = line[0].strip()
		if tempkey[:3] == tempkey[4:7]:
			tempkey = tempkey[4:]
		if tempkey not in chromids.keys():
			#ID, Sequence
			chromids[tempkey] = [[line[0], adjseq, line[4], line[3]]]
		else:
			chromids[tempkey].append([line[0], adjseq, line[4], line[3]])
	
f.close()
textfile.close()

#*******PRINT BREAK
print 'Creating ID list...'

#******************************READING INPUT FASTA FILE***************************
seqids = []

infasta = open('tempfasta.fa', 'r')
for line in infasta:
	if line[:1] == '>':
		tempid = line[1:].strip()
		if tempid[:3] == tempid[4:7]:
			tempid = tempid[4:]
		seqids.append(tempid)

infasta.close()

#*******PRINT BREAK
print 'Blasting database...'

#***************************CREATING XML FILE*******************************
#blastn -db testdb -query chunks.fa -out results.xml -outfmt 5
call(['blastn', '-db', database, '-query', 'tempfasta.fa', '-num_alignments', '8', '-out', 'temp.xml', '-outfmt', '5']) #formerly -ungapped
call(["rm", "tempfasta.fa"])

#*****************************IUPAC CHECKER METHOD**********************************
specialchars = ['R', 'Y', 'S', 'W', 'K', 'M'] #new condition
def hasSNP(sequence):
	for item in specialchars:
		if item in sequence:
			return True
	return False

#*************************SEQUENCE SNP FINDER METHOD***************************
unmap = {'R':'G/A', 'Y':'C/T', 'S':'G/C', 'W':'A/T', 'K':'G/T', 'M':'T/G'}
def snpfinder(sequence):
	for i in range(0, len(sequence)):
		if sequence[i] in specialchars:
			return [str(i), unmap[sequence[i]]]
	return None
			

#**************************V6 XML REDUCER********************************
keytags = ['<Iteration_query-def>','<Hit_id>','<Hit_len>','<Hsp_hit-from>','<Hsp_hit-to>','<Hit_def>','<Hit_num>','<Hsp_positive>','<Hsp_qseq>','<Hsp_hseq>','<Hsp_midline>','</Hit>']
infile = open('temp.xml', 'r')
outfile = open('tempsmall.xml', 'w')

print 'Reducing XML file...'
xmltotal = 0
for line in infile:
	if any(keytag in line for keytag in keytags):
		outfile.write(line)
		xmltotal += 1

infile.close()
outfile.close()

#*******PRINT BREAK
print 'Reticulating XML file...'

#************************MAIN LOOP*******************************
infile = open('tempsmall.xml', 'r')
report = open('report.csv', 'w')

commands = []
currid = ''
currmin = 0
currmax = 0
currlen = 0
legit = False
currseq = ''
idorder = []
currhit = ''
hitnum = 0
hspscore = ""
first = True
hseq = ""
midline = ""
changelocs = []
changes = []
#hitcount = 0#DEBUGLINE
#successcount = 0#DEBUG LINE
#goodcount = 0#DEBUG LINE
#weirdcount = 0#DEBUG LINE
#unfoundcount = 0#DEBUG LINE
currname = ''#V5 ADDITION
xmlcount = 0#V6 ADDITION
for line in infile:
	if '<Iteration_query-def>' in line:
		currname = line.strip().replace('<Iteration_query-def>','').replace('</Iteration_query-def>','')
	if '<Hit_id>' in line:
		currid = line.strip().replace('<Hit_id>','').replace('</Hit_id>','')
	if '<Hit_len>' in line:
		currlen = int(line.strip().replace('<Hit_len>','').replace('</Hit_len>',''))
	if '<Hsp_hit-from>' in line:
		currmin = int(line.strip().replace('<Hsp_hit-from>','').replace('</Hsp_hit-from>',''))
	if '<Hsp_hit-to>' in line:
		currmax = int(line.strip().replace('<Hsp_hit-to>','').replace('</Hsp_hit-to>',''))
	if '<Hit_def>' in line:
		for seqid in seqids:
			strippedline = line.strip().replace('<Hit_def>','').replace('</Hit_def>','')
			currhit = strippedline
			if strippedline == seqid:
				legit = True
	if '<Hit_num>' in line:
		hitnum = int(line.strip().replace('<Hit_num>','').replace('</Hit_num>',''))
	if '<Hsp_positive>' in line:
		hspscore = int(line.strip().replace('<Hsp_positive>','').replace('</Hsp_positive>',''))
	if '<Hsp_qseq>' in line:
		currseq = line.strip().replace('<Hsp_qseq>','').replace('</Hsp_qseq>','')
	#<Hsp_hseq> <Hsp_midline>
	if '<Hsp_hseq>' in line:
		hseq = line.strip().replace('<Hsp_hseq>','').replace('</Hsp_hseq>','')
	if '<Hsp_midline>' in line:
		midline = line.strip().replace('<Hsp_midline>','').replace('</Hsp_midline>','')
		#set of ' ' indices
		changelocs = []
		changes = []
		for i in range(len(midline)):
			if midline[i] == ' ':
				changelocs.append(i)
		for loc in changelocs:
			changes.append(str(loc + 1) + ':' +currseq[loc] + '/' + hseq[loc])			

	if '</Hit>' in line and hitnum < 21 and hspscore > 58: #formerly 0
		complement = '' #reverse complement testing
		if currmin > currmax:
			complement = ' (RC)'
		currfrom = min(currmin, currmax) #Accounting for reverse compliment: min > max
		currto = max(currmin, currmax)
		legit = False
		found = False
		#TODO ditch all chromids.keys business, somehow. Maybe an else... just add raw currhit to idorder? Bogus unnecessary
		#if hitnum == 1:#DEBUG LINE
		#	hitcount += 1#DEBUG LINE
		if (currname in chromids.keys()) and (hitnum == 1):
			for possible in chromids[currname]:
				if possible[1].strip() == currseq and found == False:
					templen = currto - currfrom #For printing snp location
					tempbegin = buff + min(currfrom - buff, 0) #For printing snp location
					idorder.append(possible[0] + ',' + possible[1] + ',' + possible[2] + ',' + possible[3])
					found = True
					#report filling
					if hitnum == 1 and not first:
						report.write('\n')
					first = False
					report.write(',"' + currhit[0:3] + '(' + str(hspscore) + ')' + '~' + ';'.join(changes) + '"')
					#successcount += 1#DEBUG LINE
			if found == False: #debug line TODO TODO TODO this may be plausible to salvage
				#unfoundcount += 1#debug line
				
				if snpfinder(currseq) != None: #TODO TODO TODO NEW STUFFFF doesn't woooooorkkkk. Or it does??
					try:
						idorder.append('"' + str(currname) + '",' + currseq + ',' + ','.join(snpfinder(currseq)))
					except(TypeError):
						print "I AM ERROR"
						print currhit
						print currseq
						exit(1)
					#report filling
					if hitnum == 1 and not first:
						report.write('\n')
					first = False
					bogus = False
					report.write(',"' + currhit[0:3] + '(' + str(hspscore) + ')' + '~' + ';'.join(changes) + '"')
		elif hitnum == 1:
			if hasSNP(currseq): #debug lines TODO TODO TODO
				#goodcount += 1

				idorder.append('"' + str(currhit) + '",' + currseq + ',' + ','.join(snpfinder(currseq))) #TODO TODO TODO test this
				#report filling
				if hitnum == 1 and not first:
					report.write('\n')
				first = False
				bogus = False
				report.write(',"' + currhit[0:3] + '(' + str(hspscore) + ')' + '~' + ';'.join(changes) + '"')
			#else: #TODO TODO TODO 
			#	weirdcount += 1

		elif hitnum > 1:
			report.write(',"' + currhit[0:3] + '(' + str(hspscore) + ')' + '~' + ';'.join(changes) + '"')
		#if found == True:
		#	commands.append([currid, max(currfrom - buff, 1), min(currto + buff, currlen)])
	
	#V6 additions...
	xmlcount += 1
	if xmlcount % 100 == 0:
		print 'line ' + str(xmlcount) + ' of ' + str(xmltotal)

infile.close()
report.close()

#call(["rm", "temp.xml"])

#print 'total ' + str(hitcount)#debugline
#print 'success ' + str(successcount)#debug line
#print 'unfound ' + str(unfoundcount)#debug line
#print 'good ' + str(goodcount)#debug line
#print 'weird ' + str(weirdcount)#debug line
#print 'chromids: ' + str(len(chromids))#Debug line
#print 'COMMANDS: ' + str(len(commands)) #Debug line
#print 'IDS: ' + str(len(idorder)) #Debug line

#*****************************ADDING THE IDs*******************************************
f = open('reportx.csv', 'w')
idtemp = open('report.csv', 'r')

index = 0
for line in idtemp:
	if index > len(idorder) - 1:#Debugging: this should not trigger
		print 'UNHELPFUL ERROR MESSAGE: ID MATCH FAILURE'
	else:
		f.write(idorder[index] + line)
	index += 1

idtemp.close()
f.close()
call(["rm", "report.csv"])

#*****************************IUPAC DE-MAPPER METHOD**********************************
specials = ['R', 'Y', 'S', 'W', 'K', 'M'] #new condition
demap = {'R':['G','A'], 'Y':['C','T'], 'S':['G','C'], 'W':['A','T'], 'K':['G','T'], 'M':['T','G'], 'G':['G'], 'A':['A'], 'T':['T'], 'C':['C']}
def iuUnPac(orig, new):
	if orig == new:
		return True
	if new.split('/')[0] in specials: #new condition
		if (orig.split('/')[0] in demap[new.split('/')[0]]) and (orig.split('/')[1] == new.split('/')[1] or (new.split('/')[1] in demap[new.split('/')[0]])): #crazy bool
			return True
	return False 

#*************************SAME/DIFF ANALYSIS*******************************************
f = open('reportxx.csv', 'w')
difftemp = open('reportx.csv', 'r')

csvreader = csv.reader(difftemp, delimiter=',', quotechar='"')
for line in csvreader:
	if len(line) > 4:
		for i in range(4, len(line)):
			same = 'd'
			chrompair = line[i].split('~')
			if len(chrompair) > 1 and chrompair[1] != '':
				changelist = chrompair[1].split(';')
				if len(changelist) < 5:
					for item in changelist:
						specs = item.split(':')
						#print specs #debug line
						if iuUnPac(line[3].strip(), specs[1].strip()):#(specs[1].strip() == line[3].strip())): #old location match: (specs[0].strip() == line[2].strip())
							#print specs[0] + ',' + line[2] + ',' + specs[1] + ',' + line[3] #debug line
							same = 's'

			line[i] = chrompair[0]  + ':' + same
	f.write('"' + '","'.join(line) + '"\n')

f.close()
difftemp.close()
#call(["rm", "reportx.csv"])			

#******************************THE END*****************************************
print 'FINISHED! New file is reportxx.csv'
