#Chromosome      Start   End     Probe_Number    Segment_Mean
#1       554268  639581  3       -0.561
#1       736483  9985450 1500    0.5427
#1       9987720 10304808        95      0.2286
#1       10306679        10941525        115     0.5325
#1       10946936        15407023        651     -0.3389
#1       15409548        25457903        1816    -0.5529
#1       25482036        25529549        4       -1.6003
#1       25538212        27483586        409     -0.5603

import sys

if len(sys.argv)!=5:
	print sys.argv[0],"fname length parts spike"
	sys.exit(2)

fn=sys.argv[1]
l=int(sys.argv[2])
p=int(sys.argv[3])
s=float(sys.argv[4])


segp=0

h=open(sys.argv[1])
for line in h:
	if line[0]=='C':
		continue
	line=line.strip().split()
	segl=int(line[2])-int(line[1])
	segs=float(line[4])
	if segl>l and segs>s:
		segp+=1	
h.close()

if segp<p:
	print "F",sys.argv[1]
	sys.exit(1)
else:
	print "T",sys.argv[1]
	sys.exit(0)
