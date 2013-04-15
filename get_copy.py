#!/usr/bin/python


import sys



bps=set()

its=[]

for line in sys.stdin:
	chra,coorda,chrb,coordb = map(int, line.split('\t'))
	a=(chra,coorda)
	b=(chrb,coordb)
	bps.add(a)
	bps.add(b)
	if a<b:
		its.append((a,b))
	else:
		its.append((b,a))
	
bps=list(bps)
bps.sort()

d={}

for a,b in its:
	for x in range(1,len(bps)):
		p = bps[x]
		e = (bps[x-1],bps[x])
		if a<p and p<=b:
			if e not in d:
				d[e]=0
			d[e]+=1


for x in range(1,len(bps)):
	e = (bps[x-1],bps[x])
	s="%d\t%d ~ %d\t%d\t" % (bps[x-1][0],bps[x-1][1],bps[x][1],abs(bps[x][1]-bps[x-1][1]))
	if e in d:
		print s,d[e]
	else:
		print s,0
