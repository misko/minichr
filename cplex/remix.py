import sys

simple_contigs_fname=sys.argv[1]                                                                                                                                                                                                      
                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                      
simple_contigs=[]                                                                                                                                                                                                                     
simple_contigs_f=open(simple_contigs_fname)                                                                                                                                                                                           
for line in simple_contigs_f:                                                                                                                                                                                                         
        simple_contigs.append(eval(line.strip())[0])

flows={}

for line in sys.stdin:
	if line[0]=='c':
		line=line.strip().split()
		i=int(line[0][1:])
		m=int(line[1])
		for x in range(len(simple_contigs[i])):
			fn=simple_contigs[i][x]
			tn=simple_contigs[i][(x+1)%len(simple_contigs[i])]
			if fn not in flows:
				flows[fn]={}
			if tn not in flows[fn]:
				flows[fn][tn]=0
			flows[fn][tn]+=m

for fn in flows:
	for tn in flows[fn]:
		print "f\t%d\t%d\t%d" % (fn,tn,flows[fn][tn])
		
		
