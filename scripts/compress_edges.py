#!/usr/bin/bash


import sys


arcs={}
nodes=[]


comment=""
for line in sys.stdin:
	if line[0]=='c':
		comment=line.strip()
	elif line[0]=='a':
		x=tuple(map(int,line.split()[1:3]))
		if x not in arcs:
			arcs[x]=[]
		if len(comment)!=0:
			arcs[x].append(comment)
		comment=""
		arcs[x].append(line.strip())
	elif line[0]=='n':
		if len(comment)!=0:
			nodes.append(comment)
		comment=""
		nodes.append(line.strip())
