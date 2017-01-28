#!/usr/bin/env python

import re
import sys
from collections import defaultdict

old="/Users/kristen/Documents/transposon_figure_data/data/compare_maps_old.txt"
new="/Users/kristen/Documents/transposon_figure_data/data/compare_maps_new.txt"


comp=defaultdict(list)
with open(old, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)

		trait=items[0]
		for index,value in enumerate(items[1:153]):
			tt=(str(index),str(value))
			tog=":".join(tt)
			comp[trait].append(tog)

OUT=open("/Users/kristen/Documents/transposon_figure_data/data/compare_differences.txt", 'w')
OUT1=open("/Users/kristen/Documents/transposon_figure_data/data/compare_values.txt", 'w')
with open(new, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait=items[0]
		old_dat=comp[trait]
		OUT.write(trait)
		OUT1.write(trait)
		for index,value in enumerate(items[1:153]):
			old_info=old_dat[index]
			old_items=re.split(":",old_info)
			old_index=old_items[0]
			old_value=old_items[1]
			value_diff=int(value)-int(old_value)
			OUT.write('\t' + str(value_diff))
			OUT1.write('\t' +str(old_value) + "/" + str(value))
		OUT.write('\n')
		OUT1.write('\n')

OUT.close()
OUT1.close()