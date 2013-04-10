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

def from_chr(x):
	if x==23:
		return 'chrX'
	if x==24:
		return 'chrY'
	if x==25:
		return 'chrM'
	return "chr%d" % x

from math import sqrt

h=open(arc_coverage_filename,'r')
for line in h:
	line=line.split()
	chr=to_chr(line[0])
	coord=int(line[1])
	support=int(line[2])	
	arcs[(chr,coord)]=float(support)-sqrt(float(support))


for line in sys.stdin:
	if line[:3]!="chr":
		continue
	line=line.split()
	#chr1:59446789+  chr1:59446926+  13
	
	chra=to_chr(line[0].split(':')[0])
	coorda=int(line[0].split(':')[1].replace('+','').replace('-',''))
	astrand='+'
	if '-' in line[0]:
		astrand='-'
	a=(chra,coorda)	

	chrb=to_chr(line[1].split(':')[0])
	coordb=int(line[1].split(':')[1].replace('+','').replace('-',''))
	bstrand='+'
	if '-' in line[1]:
		bstrand='-'
	b=(chrb,coordb)

	s=int(line[2])
	#s=s+sqrt(s)

	if a in arcs and b in arcs:
		print "%s\t%d\t%s\t%s\t%d\t%s\t%d" % (from_chr(chra),coorda,astrand,from_chr(chrb),coordb,bstrand,round(s/((arcs[a]+arcs[b])/4)))


