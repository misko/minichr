#!/usr/bin/python


import sys

if len(sys.argv)!=2:
	print "%s collapse_size" % sys.argv[0]
	sys.exit(1)

collapse_size=int(sys.argv[1])

d={}

for line in sys.stdin:
	s=line.split()
	#chr1    100992522       47591589        0       2       2       31      4.39049e+08     6873    chr3    EDGE
	#chr1    101283258       101283591       0       2       6       3       438     6888    chr1    EDGE
	#chr1    101719157       83133391        1       2       52      1       -1.43252e+09    2489    chrX    EDGE
	#chr1    101724032       101724287       0       3       78      90      436.5   6909    chr1    EDGE
	#chr1    102018665       102018985       0       2       13      26      439.5   6918    chr1    EDGE
	chra=s[0]
	posa=int(int(s[1])/collapse_size)*collapse_size
	chrb=s[-2]
	posb=int(int(s[2])/collapse_size)*collapse_size
	p1=(chra,posa)
	p2=(chrb,posb)
	ty=int(s[3])
	if p1>p2:
		t=p1
		p1=p2
		p2=t
	k=(p1,p2,ty)
	if k not in d:
		d[k]=0
	d[k]+=int(s[4])

for k in d:
	p1,p2,ty=k
	chra,posa=p1
	chrb,posb=p2
	print "%s\t%d\t%d\t%d" % (chra,posa,posb,ty),
	print "%d\t0\t0\t0.0\t0\t%s\tEDGE" % (d[k],chrb)


