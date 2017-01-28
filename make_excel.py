#!/usr/bin/env python
import sys
import re
import codecs
import xlsxwriter
import os

infile=sys.argv[1]
out_name=os.path.splitext(os.path.basename(infile))[0]

# Open the input file with the correct encoding.
textfile = codecs.open(infile, 'r', 'utf-8')

# Create an new Excel file and convert the text data.
workbook = xlsxwriter.Workbook(out_name + '.xlsx')
worksheet = workbook.add_worksheet()

# Widen columns to make the text clearer.
worksheet.set_column('A:A', 25)
worksheet.set_column('B:B', 20)
worksheet.set_column('C:C', 20)

bold = workbook.add_format({'bold': True})
# Start from the first cell.
row = 0
col = 0

# Read the text file and write it to the worksheet.
for line in textfile:
	line=line.rstrip('\n')
	col=0
	items=re.split('\t',line)
	for i in items:
		if row==0:
			worksheet.write(row, col,i,bold)
		else:
			worksheet.write(row, col,i)
		col+=1
	row += 1

workbook.close()
textfile.close()