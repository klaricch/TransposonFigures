#!/usr/bin/env python
import re
from subprocess import Popen, PIPE

file1="/Users/kristen/Documents/transposon_figure_data/tables/TajimaD_Table.txt"
file2="/Users/kristen/Documents/transposon_figure_data/tables/taj_genes.txt"

OUT=open("/Users/kristen/Documents/transposon_figure_data/figures/tmp.txt",'w')

genes={}
with open(file2, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		ID,gene=items[0:2]
		genes[ID]=gene

with open(file1, 'r') as IN:
	headers=next(IN).rstrip('\n')
	headers=headers+'\t'+'Genes in Bin' + '\n'
	OUT.write(headers)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		Chrom,Bin_Start,Number_TE,TJ=items[0:4]
		ID=Chrom+Bin_Start
		ID_genes=genes[ID]
		OUT.write("{Chrom}\t{Bin_Start}\t{Number_TE}\t{TJ}\t{ID_genes}\n".format(**locals()))

OUT.close()

result, err = Popen(["""mv /Users/kristen/Documents/transposon_figure_data/figures/tmp.txt /Users/kristen/Documents/transposon_figure_data/tables/TajimaD_Table.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()		
