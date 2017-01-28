#!/usr/bin/env python
#this script fixes TE names in count matrix

import re
from subprocess import Popen, PIPE
file1="/Users/kristen/Documents/transposon_figure_data/tables/T_kin_C_matrix_full_reduced.txt"
OUT=open("/Users/kristen/Documents/transposon_figure_data/figures/tmp.txt",'w')
with open(file1) as IN:
	headers=next(IN)
	OUT.write(headers)
	for line in IN:
		if re.search("^coverage",line):
			OUT.write(line)
		else:
			line=line.rstrip('\n')

			items=re.split('\t',line)
			te_info=items[0]
			line_info='\t'.join(items[1:])
			
			more_info=re.split('_TRANS_',te_info)
			method,te=more_info[0:2]
			te=re.sub("_C$","",te) #strip trailing _C
			if method!="reference" and method!="ZERO_new":
				if method=="cumulative":
					abbrev="(cu)"
				elif method=="ONE_new" or method=="new":
					abbrev="(ins)"
				elif method=="absent":
					abbrev="(AR)"
				new_name=te+abbrev
				OUT.write("{new_name}\t{line_info}\n".format(**locals()))
OUT.close()
result, err = Popen(["""mv /Users/kristen/Documents/transposon_figure_data/figures/tmp.txt /Users/kristen/Documents/transposon_figure_data/tables/tmp_T_kin_C_matrix_full_reduced.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()		


