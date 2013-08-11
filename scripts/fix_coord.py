#!/usr/bin/python



import sys

def getstart(s):
	c=s.split(":")[0]
	s=s.split(":")[1].split("-")[0]
	return (c,int(s))

for line in sys.stdin:
	if line[0]=="@":
		print line,
		continue
	line=line.split()
	if line[2]!="*":
		s1=getstart(line[2])
		line[2]=s1[0]
		line[3]=str(int(line[3])+s1[1])
		if line[6]!="*":
			s2=s1
			if line[6]!='=':
				s2=getstart(line[6])
			line[6]=s2[0]
			line[7]=str(int(line[7])+s2[1])
		line[8]="77777"
	print "\t".join(line)
			
				
