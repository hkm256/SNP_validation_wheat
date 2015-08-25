import csv

infile = open('reportxx.csv', 'r')
goodfile = open('reportxx_good.csv', 'w')
badfile = open('reportxx_bad.csv', 'w')
csvreader = csv.reader(infile, delimiter=',', quotechar='"')

for row in csvreader:
	if len(row) > 4:
		chroms = []
		for i in range(4, len(row)):
			chunks = row[i].split('(')
			if 's' in chunks[1] and chunks[1][0] == '6':
				chroms.append(chunks[0])
		chroms = set(chroms)
		if len(chroms) > 1:
			badfile.write('"' + '","'.join(row) + '"\n')
		else:
			goodfile.write('"' + '","'.join(row) + '"\n')

infile.close()
goodfile.close()
badfile.close()
