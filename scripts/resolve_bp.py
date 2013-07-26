#!/usr/bin/python



import sys


if len(sys.argv)!=4:
	print "%s assembled.fa ref1.fa ref2.fa" % sys.argv[0]
	sys.exit(1)


def reverse_complement(s):
	s2=""
	for x in s:
		if x=="c":
			s2="g"+s2
		elif x=="C":
			s2="G"+s2
		elif x=="a":
			s2="t"+s2
		elif x=="A":
			s2="T"+s2
		elif x=="g":
			s2="c"+s2
		elif x=="G":
			s2="C"+s2
		elif x=="t":
			s2="a"+s2
		elif x=="T":
			s2="A"+s2
		elif x=="N":
			s2="N"+s2
		else:
			print "ERROR"
			sys.exit(1)
	return s2


def read_fa(f):
	h=open(f,'r')
	lines=h.readlines()
	h.close()
	return (lines[0][1:].strip(),"".join(map(lambda x : x.strip(),lines[1:])))

afa=read_fa(sys.argv[1])
r1fa=read_fa(sys.argv[2])
r2fa=read_fa(sys.argv[3])

r1_algn=()
r2_algn=()


mxa=0
mxb=0
for alignment in sys.stdin.readlines():
	alignment=alignment.strip().split()
	sc=int(alignment[0])
	if sc<150:# or int(alignment[0])<=mx:
		#print >> sys.stderr, "#SKIPPING"
		continue
	#it blat
	#956     1       0       0       0       0       0       0       +       Genome_LEN_2324 2324    1306    2263    chr7:54413327-54414627  1301    344     1301    1       957,    1306,   344,
	#987     0       0       0       0       0       0       0       +       Genome_LEN_2324 2324    320     1307    chr12:60845840-60847140 1301    0       987     1       987,    320,    0,
	a=[alignment[9],int(alignment[11])+1,int(alignment[12])+1]
	b=[alignment[13],int(alignment[15])+1,int(alignment[16])+1]
	b_reverse=False
	if alignment[8]=='-':
		b_reverse=True
		t=b[1]
		b[1]=b[2]
		b[2]=t
		

	if b[0]==r1fa[0]:
		if sc>mxa:
			mxa=sc
			r1_algn=(a,b,b_reverse)
	elif b[0]==r2fa[0]:
		if sc>mxb:
			mxb=sc
			r2_algn=(a,b,b_reverse)
	else:
		print >> sys.stderr, "ERROR!"
	#print b

if len(r1_algn)==0 or len(r2_algn)==0:
	print "#-1\t#no bp"
	sys.exit(1)

#print r1_algn, r2_algn

if r2_algn[0][1]<r1_algn[0][1]:
	t=r2_algn
	r2_algn=r1_algn
	r1_algn=t
	t=r1fa
	r1fa=r2fa
	r2fa=t

insert=[]
micro_homology=[]
if r1_algn[0][2]>r2_algn[0][1]:
	sz=r1_algn[0][2]-r2_algn[0][1]
	if not r1_algn[2]:
		micro_homology.append(r1fa[1][r1_algn[1][1]-1:r1_algn[1][2]-1][-sz:])
	else:
		micro_homology.append(reverse_complement(r1fa[1][r1_algn[1][2]-1:r1_algn[1][1]-1])[-sz:])
	if not r2_algn[2]:
		micro_homology.append(r2fa[1][r2_algn[1][1]-1:r2_algn[1][2]-1][:sz])
	else:
		micro_homology.append(reverse_complement(r2fa[1][r2_algn[1][2]-1:r2_algn[1][1]-1])[:sz])
elif r1_algn[0][2]<r2_algn[0][1]:
	insert.append(afa[r1_algn[0][2]-1:r2_algn[0][1]-1])


def to_pos(s):
	chr=s.split(':')[0]
	f,t=map(lambda x : int(x) ,s.split(':')[1].split('-'))
	return [chr,f,t+1]

r1p=to_pos(r1fa[0])
r2p=to_pos(r2fa[0])

#r1 is first, followed by r2
if not r1_algn[2] and not r2_algn[2]:
	# both forward
	#type 0
	ty=0
	#print r1p,r1_algn[1]
	r1p.append(r1p[1]+r1_algn[1][2]-1)
	r2p.append(r2p[1]+r2_algn[1][1]-1)
elif r1_algn[2] and not r2_algn[2]:
	#goes from back to forward
	# type 3
	r1p.append(r1p[1]+r1_algn[1][2]-1)
	r2p.append(r2p[1]+r2_algn[1][1]-1)
	ty=3
elif not r1_algn[2] and r2_algn[2]:
	#goes from forward to back
	r1p.append(r1p[1]+r1_algn[1][2]-1)
	r2p.append(r2p[1]+r2_algn[1][1]-1)
	# type 2
	ty=2
elif r1_algn[2] and r2_algn[2]:
	#goes from back to back
	#type 1
	r1p.append(r1p[1]+r1_algn[1][2]-1)
	r2p.append(r2p[1]+r2_algn[1][1]-1)
	ty=1


print r1p[0],r1p[3],r2p[0],r2p[3],ty,insert,micro_homology

#	print a[1],a[2],afa[1][a[1]-1:a[2]]
#	if b_reverse:
#		print b[1],b[2],reverse_complement(ref[1][b[2]-1:b[1]])
#	else:
#		print b[1],b[2],ref[1][b[1]-1:b[2]]
#	print a,b
	

#Genome_LEN_1793 1277 1563; chr12:60845840-60847140 987 701; score = 27001.000000 (-)
#Genome_LEN_1793 321 1280; chr7:54413327-54414627 1301 342; score = 90161.000000 (-)

