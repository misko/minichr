#!/usr/bin/python


import sys
import gzip

if len(sys.argv)!=2:
	print "%s fname_to_remove" % sys.argv[0]
	sys.exit(1)

rm=set()

h=gzip.open(sys.argv[1],'r')
for line in h:
	line=line.strip()
	rm.add(line)

for line in sys.stdin:
	line=line.strip()
	if line not in rm:
		print line
