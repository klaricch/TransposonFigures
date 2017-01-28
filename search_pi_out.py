#!/usr/bin/evv python
# this script check is any of the variation in pirnas in strains with outlier te familes is in the pirna aligning to that te
import re
from collections import defaultdict
import sys
pi_matches=defaultdict(list)
with open("/Users/kristen/Documents/transposon_figure_data/tables/table_pirnas.txt",'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		transcript,te=items[0:2]
		method=items[3]
		if method=="Both" or method=="BLAST Only":
			pi_matches[transcript].append(te)


#sys.exit()
with open("/Users/kristen/Documents/transposon_figure_data/figures/pi_outliers.txt",'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		transcript=items[4]
		tes=items[7]
		te_fams=tes.split(",")
		te_fams=[re.sub("\(.*\)","",i) for i in te_fams]

		if transcript in pi_matches:

			compare=pi_matches[transcript]
			#print transcript
			#print compare
			for i in te_fams:
				#print i
				
				if i in compare:
				#if re.search(i,compare):
					print "FOUND"
					print i
