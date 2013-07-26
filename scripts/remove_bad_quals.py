#!/usr/bin/python

import sys


if len(sys.argv)!=5:
	print "%s f1fq f2fq o1 o2" % sys.argv[0]


def getreads(f):
	h=open(f,'r')
	reads=[]
	read=[]
	print "reading reads"
	for line in h:
		read.append(line.strip())
		if len(read)==4:
			reads.append(read)
			read=[]
	return reads	



fr=getreads(sys.argv[1])
rr=getreads(sys.argv[2])

o1=open(sys.argv[3],'a')
o2=open(sys.argv[4],'a')

for x in range(len(fr)):
	rfr=fr[x]
	rrr=rr[x]
	if len(rfr[1].replace('n','').replace('N',''))<30:
		pass
	elif len(rrr[1].replace('n','').replace('N',''))<30:
		pass
	else:
		o1.write("\n".join(rfr)+"\n")
		o2.write("\n".join(rrr)+"\n")

