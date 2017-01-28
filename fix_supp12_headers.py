#!/usr/bin/env python

import re
from subprocess import Popen, PIPE
file1="/Users/kristen/Documents/transposon_figure_data/tables/T_with_monomorphic_kin_matrix_full.txt"
OUT=open("/Users/kristen/Documents/transposon_figure_data/figures/tmp.txt",'w')
with open(file1) as IN:
	header=next(IN)
	header=header.rstrip()
	for i in re.split('\t',header):
		print "start:"
		print i
		if i == "Marker":
			OUT.write("Strain")
		elif re.search("_non-reference_",i):
			i=re.sub("_non-reference_(.*?)$","_NR",i)
			OUT.write('\t' + i)
		else:
			i=re.sub("$","_R\t",i)
			print i

			if len(i)>3:
				print "YES"
				OUT.write('\t' + i)
		print 'hold'

	print "end"
	OUT.write('\n')
	for line in IN:
		OUT.write(line)
OUT.close()
result, err = Popen(["""mv /Users/kristen/Documents/transposon_figure_data/figures/tmp.txt /Users/kristen/Documents/transposon_figure_data/tables/tmp_T_with_monomorphic_kin_matrix_full.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()		


#file2="/Users/kristen/Documents/transposon_figure_data/tables/T_kin_matrix_full.txt"
#OUT=open("/Users/kristen/Documents/transposon_figure_data/figures/tmp.txt",'w')
#with open(file2) as IN:
#	header=next(IN)
#	for i in re.split('\t',header):
#		if i == "Marker":
#			OUT.write("Strain")
#		elif re.search("_non-reference_",i):
#			i=re.sub("_non-reference_(.*?)$","_NR",i)
#			OUT.write('\t' + i)
#		else:
#			i=re.sub("$","_R\t",i)
#			OUT.write('\t' + i)
#	OUT.write('\n')
#	for line in IN:
#		OUT.write(line)
#OUT.close()
#result, err = Popen(["""mv /Users/kristen/Documents/transposon_figure_data/figures/tmp.txt /Users/kristen/Documents/transposon_figure_data/tables/tmp_T_kin_matrix_full.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()		
