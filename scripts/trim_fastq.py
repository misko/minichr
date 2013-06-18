#!/usr/bin/python

import sys


h=open(sys.argv[1],'r')



lines=[]
for line in h:
	lines.append(line.strip())
	if len(lines)==4:
		i=0
		while i<len(lines[3]) and lines[3][len(lines[3])-i-1]=="#":
			i=i+1
		if len(lines[3])-i>10:
			print lines[0]
			print lines[1][:-i-1]
			print lines[2]
			print lines[3][:-i-1]
		else:
			print lines[0]
			#print lines[1][:-i-1]
			print "N"*80
			print lines[2]
			#print lines[3][:-i-1]
			print "#"*80
		lines=[]	

