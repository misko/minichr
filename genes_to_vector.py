#!/usr/bin/python

import sys


if len(sys.argv)!=2:
	print "%s gene list"  % sys.argv[0]
	sys.exit(1)

glf=sys.argv[1]

gl=set()

h=open(glf)
for line in h.readlines():
	line=line.strip().replace('*','')
	if len(line)>1:
		gl.add(line)
h.close()

gl=list(gl)
gl.sort()

patients={}
labels={}
lk=set()
bins=[72,30,15,6]
for line in sys.stdin:
	line=line.strip().split()
	id=line[0]
	label=line[1]
	lk.add(label)
	if id not in patients:
		patients[id]={}
	if id not in labels:
		labels[id]=label
	for x in range(4):
		l=map(lambda x : x.replace('*',''), line[x+3].split(','))
		for g in l:
			if g in gl:
				if g not in patients[id]:
					patients[id][g]=bins[x]
lk=list(lk)
lk.sort()
pk=patients.keys()
pk.sort()
s=[]
s.append('id')
s.append('label')
for g in gl:
	s.append(g)
#print ",".join(s)

for p in patients:
	v=[str(pk.index(p)),str(lk.index(labels[p]))]
	for x in range(len(gl)):
		if gl[x] in patients[p]:
			v.append(str(patients[p][gl[x]]))
		else:
			v.append("0")
	print ",".join(v)

h=open('lookup','w')
print >> h, lk
print >> h, pk
h.close()
	
