#!/usr/bin/evv python
import re
import sys
close_dict={}
with open("/Users/kristen/Documents/transposon_figure_data/tables/close_table.txt",'r') as IN:
	next(IN)
	for line in IN:

		line=line.rstrip()
		items=line.split('\t')
		te_id=items[0]+"_"+items[1]
		if not re.search("(ref)",te_id):
			close_dict[te_id]=0
print len(close_dict)



median_dict={}
with open("/Users/kristen/Documents/transposon_figure_data/tables/median_table.txt",'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te_id=items[0]+"_"+items[1]
		if not re.search("(ref)",te_id):
			median_dict[te_id]=0
print len(median_dict)
#sys.exit()
counta=0
countb=0

OUT=open("tmp.txt", 'w')
with open("/Users/kristen/Documents/transposon_figure_data/tables/Peak_Table.txt",'r') as IN:
	header=next(IN)
	header=re.sub('\n','\t',header)
	OUT.write("{header}Distance to TE Position\tMedian Value\n".format(**locals()))
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te_id=items[0]+"_"+items[1]
		if te_id in close_dict:
			counta+=1
			distance="near"
		elif re.search("total",te_id):
			distance="NA"
		else:
			distance="far"

		if te_id in median_dict:
			countb+=1
			median="same"
		else:
			median="different"
		OUT.write("{line}\t{distance}\t{median}\n".format(**locals()))

print counta
print countb
OUT.close()