#/usr/bin/python

import gzip
import sys
from bisect import bisect_left

if len(sys.argv)<4:
	print "%s n master f1 f2 ..." % sys.argv[0]
	sys.exit(1)

n=int(sys.argv[1])
master_filename=sys.argv[2]
input_filenames=sys.argv[3:]

h=[0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]

for i in range(1,len(h)):
	h[i]+=h[i-1]

def chr_to_int(s):
	if s[:3]=="chr":
		s=s[3:]
	if s=="x":
		s="23"
	elif s=="y":
		s="24"
	elif s=="m":
		s="25"
	return int(s)

def to_hg(p):
	chr,spos,epos=p
	offset=h[chr-1]
	return (-1,spos+offset,epos+offset)

def line_to_pos(l):
	line=l.lower().strip().split()
	if line[2] in ('-','+') and line[5] in ('-','+'):
		#ctx format
		chr=chr_to_int(line[0])
		spos=int(line[1])
		epos=int(line[1])
		k=(chr,spos,epos)
		chr=chr_to_int(line[3])
		spos=int(line[4])
		epos=int(line[4])
		kp=(chr,spos,epos)
		return (k,kp)
	else:
		chr=chr_to_int(line[0])
		spos=int(line[1])
		epos=int(line[2])
		#chr,spos,epos=to_hg((chr,spos,epos))
		return (chr,spos,epos)

doubles={}
d=[]
master_file=gzip.open(master_filename,'r')
for line in master_file:
	if line[0]=='@' or line[0]=='#':
		continue
	z=line_to_pos(line)
	if len(z)==2:
		k,kp=z
		if not k in doubles:
			doubles[k]=set()
		doubles[k].add(kp)
		if not kp in doubles:
			doubles[kp]=set()
		doubles[kp].add(k)
		chr,spos,epos=k
		d.append((chr,epos,spos))
		chr,spos,epos=kp
		d.append((chr,epos,spos))
	else:
		chr,spos,epos=z
		d.append((chr,epos,spos))

sorted=True
for i in xrange(len(d)-1):
	if d[i]>d[i+1]:
		sorted=False
		break

def overlap(mi,i):
	mchr,mepos,mspos=mi
	chr,spos,epos=i
	if chr!=mchr:
		return 0
	d=min(epos,mepos)-max(spos,mspos)
	if d<0:
		return 0
	else:
		return d


if not sorted:
	print "sorting.."
	d.sort()
	print "done sorting"
else:
	print "already sorted"

for filename in input_filenames:
	file=gzip.open(filename,'r')
	for linex in file:
		if linex[0]=='@' or linex[0]=='#':
			continue
		#line=linex.lower().strip().split()
		#chr=chr_to_int(line[0])
		#spos=int(line[1])
		#epos=int(line[2])
		#chr,spos,epos=to_hg((chr,spos,epos))
		#k=(chr,spos,epos)
		z=line_to_pos(linex)
		if len(z)!=2:
			k=z
			chr,spos,epos=k
			i=bisect_left(d,(chr,spos-n,epos))
			dont_print=False
			if i<0:
				i+=1
			if i<len(d) and d[i][0]!=chr:
				i+=1
			while i<len(d) and d[i][0]==chr and d[i][2]<n+epos: #start of m-int is at most n after end pos
				#print i,overlap(d[i],k),k,d[i]
				if overlap(d[i],k)+n>epos-spos:
					#print "IN"
					dont_print=True
					break
				i+=1
			if not dont_print:
				print linex,
		else:
			k,kp=z
			chr,spos,epos=k
			i=bisect_left(d,(chr,spos-n,epos))
			found=False
			if i<0:
				i+=1
			if i<len(d) and d[i][0]!=chr:
				i+=1
			while not found and i<len(d) and d[i][0]==chr and d[i][2]<n+epos: #start of m-int is at most n after end pos
				if d[i][2]-n<spos<d[i][2]+n:
					#need to check the other feet
					for pchr,pspos,pepos in doubles[(d[i][0],d[i][2],d[i][1])]:
						if pchr==kp[0] and pspos-n<kp[1]<pspos+n:
							found=True
							break
				i+=1
			if not found:
				print linex,



