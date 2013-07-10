#!/usr/bin/python

import sys


if len(sys.argv)!=6:
	print "%s r1 r2 or1 or2 len" % sys.argv[0]
	sys.exit(1)


r1f=sys.argv[1]
r2f=sys.argv[2]
or1f=sys.argv[3]
or2f=sys.argv[4]
l=int(sys.argv[5])


h=open(r1f,'r')
r1=map(lambda x : x.strip(), h.readlines())
h.close()


h=open(r2f,'r')
r2=map(lambda x : x.strip(), h.readlines())
h.close()

if len(r1)!=len(r2):
	print >> sys.stderr, "failed precondition!"

or1=open(or1f,'w')
or2=open(or2f,'w')

for x in range(len(r1)):
	l1=len(r1[x].replace('n','').replace('N',''))
	l2=len(r2[x].replace('n','').replace('N',''))
	if l1<l or l2<l:
		#drop this
		pass
	else:
		or1.write(r1[x]+"\n")
		or2.write(r2[x]+"\n")
