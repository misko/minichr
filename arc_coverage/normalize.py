#!/usr/bin/python



import sys


if len(sys.argv)!=2:
	print "%s arc_coverage" % sys.argv[0]
	sys.exit(1)


arc_coverage_filename=sys.argv[1]


arcs={}



def to_chr(s):
	s=s.lower()
	if s[:3]=='chr':
		s=s[3:]
	if s=='x':
		return 23
	if s=='y':
		return 24
	if s=='m':
		return 25
	return int(s)

h=open(arc_coverage_filename,'r')
for line in h:
	line=line.split()
	chr=to_chr(line[0])
	coord=int(line[1])
	support=int(line[2])	
	arcs[(chr,coord)]=support
