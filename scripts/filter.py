#!/usr/bin/python


import sys

if len(sys.argv)!=5:
	print "%s mean stddev rname out_unsorted" % sys.argv[0]
	sys.exit(1)

mean=int(float(sys.argv[1]))
stddev=int(float(sys.argv[2]))
rname=sys.argv[3]
out_unsorted=sys.argv[4]

h=open(out_unsorted,'w')

l={}

header=True
i=0
for line_f in sys.stdin:
	if header and line_f[0]=='@':
		print line_f,
		continue
	if header:
		header=False
	line=line_f.split('\t')
	flags=int(line[1])
	proper=flags&2
	isize=abs(int(line[8]))
	if isize>(mean+3*stddev):
		print >> h,rname+line_f,
	elif line[6]!=line[2] and line[6]!="=":
		print >> h,rname+line_f,
	elif proper==0:
		#this means it is not a proper pair
		#check if its totally unmapped
		#0xC = 0x4+0x8
		if flags&0xC==0xC:
			#totally unmapped, dont print it!
			pass
		else:
			name=line[0]
			if not name in l:
				l[name]=(i,line_f)
			else:
				z,m_line_f=l[name]
				print rname+m_line_f,
				print rname+line_f,
				l.pop(name)
	else:
		#proper pair
		#check the cigar and mate cigar
		indel=(line[5]!="100M")
		name=line[0]
		if not name in l:
			#dont have its pair yet, store in mem
			l[name]=(i,indel,line_f)
		else:
			#have its pair
			z,m_indel,m_line_f=l[name]
			if m_indel or indel or isize>(mean+3*stddev):
				print rname+m_line_f,
				print rname+line_f,
				pass
			l.pop(name)	
	i+=1
	if i%10000==0:
		#keys=l.keys()
		#x=[]
		#for key in keys:
		#	x.append(l[key])
		#x.sort()
		#print x[:10]
		#print len(l)
		pass
h.close()
