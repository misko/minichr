#!/usr/bin/python


import sys


seq=""
for line in sys.stdin:
	if line[0]=='>':
		print line,
	else:
		seq+=line.strip()
i=0
while i<len(seq):
	print seq[i:i+80]
	i+=80
		
