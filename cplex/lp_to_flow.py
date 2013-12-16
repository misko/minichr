import sys
import gzip


if len(sys.argv)!=3:
	print "%s edge_lookup solution" % sys.argv[0]
	sys.exit(1)

el_fname=sys.argv[1]
sol_fname=sys.argv[2]

edge_lookup=[]

h=gzip.open(el_fname)
for line in h:
	line=line.strip()
	line=line[1:-1].split(',')
	e=int(line[0])
	f=int(line[1])
	t=int(line[2])
	ty=int(line[3])
	if len(edge_lookup)!=e:
		print >> sys.stderr, " Failed!"
		sys.exit(1)
	edge_lookup.append((f,t,ty))
h.close()


flows={}

h=open(sol_fname)
for line in h:
	if line[0]=='f':
		line=line[1:-1].split()
		fl=float(line[1])
		e=int(line[0])
		fn,tn,ty=edge_lookup[e]
		if ty==0:
			pass
		elif ty==1:
			fn-=1
			tn-=1
		elif ty==2:
			tn-=1
		elif ty==3:
			fn-=1
		if fn not in flows:
			flows[fn]={}
		if tn not in flows[fn]:
			flows[fn][tn]=0
		flows[fn][tn]+=fl
h.close()

for fn in flows:
	for tn in flows[fn]:
		print "f\t%d\t%d\t%d" % (fn,tn,flows[fn][tn])
