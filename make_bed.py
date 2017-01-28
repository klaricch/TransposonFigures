#!/usr/bin/env python
# this script generates the bed file to be loaded to cendr
# NOTE: run in tables dir
import re
from subprocess import Popen, PIPE

kin="/Users/kristen/Documents/transposon_figure_data/tables/final_tables/Supplemental_Table_S1.txt"
result, err = Popen(["""bash ../R_scripts/transpose_matrix.sh {kin} working_file.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
OUT=open("tes_cender.bed",'w')


strain_list={}
with open("working_file.txt", 'r') as IN:
	header=next(IN)
	print header
	header=header.rstrip()
	header=header.split('\t')
	for i,j in enumerate(header):
		strain_list[i]=j

	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		info=items[0]
		info2=info.split("_")
		
		
		chromosome=info2[0]
		position=info2[1]
		te=info2[2:]
		print te
		transposon="_".join(te)
		
		transposon=re.sub("_R$","",transposon)
		transposon=re.sub("_NR$","",transposon)
		print transposon
		position2=int(position)+1
		for i,j in enumerate(items):
			if i !="0":
				if j == "1":
					strain=strain_list[i]
					OUT.write("{chromosome}\t{position}\t{position2}\t{strain}\t{transposon}\t+\n".format(**locals()))

OUT.close()