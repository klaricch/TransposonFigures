#!/usr/bin/env python
import re
import sys
from collections import defaultdict
kin_mono="/Users/kristen/Documents/transposon_figure_data/tables/UT_with_monomorphic_kin_matrix_full.txt"
kin="/Users/kristen/Documents/transposon_figure_data/data/kin_matrix_full.txt"
ctcp="/Users/kristen/Documents/transposon_figure_data/data/CtCp_all_nonredundant.txt"
solo="/Users/kristen/Documents/transposon_figure_data/data/solo_class.txt"

classes={}
totals=defaultdict(int)
totals_class=defaultdict(int)

#need to add in non-monomospohic
with open(solo, 'r') as IN:
	header=next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te=items[0]
		if te=="TC2":
			te="Tc2"

		classification=items[1]

		classes[te]=classification

#need to add in non-monomospohic
with open(ctcp, 'r') as IN:
	header=next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te_info=items[3]
 		match=re.search("(\w+)_(\d+)_(.*)_(non-)?reference",te_info)
		te=match.group(3)
		classification=items[7]
		te=items[0]
		classification=items[1]
		classes[te]=classification



with open(kin_mono, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te=items[0]
		strains=items[1:]
		info=te.split("_")
		chrom,pos=info[0:2]
		site=info[-1]
		transposon=info[2:-1]
		transposon="_".join(transposon)

		if transposon not in classes:
			te_class="unknown"
		else:
			te_class=classes[transposon]

		if "1" in strains: # a double check
			totals[site]+=1
			totals_class[te_class]+=1


non_mono=0
with open(kin, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		te=items[0]
		strains=items[1:]
		info=te.split("_")
		chrom,pos=info[0:2]
		site=info[-1]
		transposon=info[2:-1]
		transposon="_".join(transposon)


		if "1" in strains and site=="R": # a double check
			non_mono+=1


print non_mono
print totals["R"]
mono_sites=totals["R"]-non_mono
print mono_sites

OUT=open("Table1.txt",'w')
OUT.write("Type\tTotal\n")
overall_total=sum(totals.values())
OUT.write("Total Transposon Sites\t"+str(overall_total)+"\n")
OUT.write("Insertion Sites\t"+str(totals["NR"])+"\n")
OUT.write("Active Reference Sites\t"+str(non_mono)+"\n")
OUT.write("Monomorphic Reference Sites\t"+str(mono_sites)+"\n")
OUT.write("DNA Elements\t"+str(totals_class["dnatransposon"])+"\n")
OUT.write("Retrotransposon Elements\t"+str(totals_class["retrotransposon"])+"\n")
OUT.write("Unknown Transposon Elements\t"+str(totals_class["unknown"])+"\n")
OUT.close()